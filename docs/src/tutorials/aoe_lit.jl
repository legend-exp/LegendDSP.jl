# # A/E Pulse-Shape Parameter

# Single-site events (SSE) deposit their energy in a small volume of the detector. Their current pulse is short and
# its maximum A is proportional to the energy E. Multi-site events (MSE), in which the energy is deposited
# at several places, produce a sum of pulses that arrive at different times. The maximum of the summed
# current is lower for the same energy. The ratio A/E therefore separates SSE from MSE.
# This tutorial computes A/E from simulated waveforms and compares alternative current filters.
# The charge-drift parameter Q-drift and the late-charge parameter LQ are computed as well.

# The tutorial continues the [Energy Reconstruction](@ref) tutorial and uses the same steps for
# the baseline, the pole-zero correction and the energy, without repeating their explanation.
# The calibration of A/E and the cut on it are not part of LegendDSP.jl. They are done with
# [LegendSpecFits.jl](https://github.com/legend-exp/LegendSpecFits.jl); see the
# [documentation](https://legend-exp.github.io/LegendSpecFits.jl/stable/) of that package.

# ## Packages

using ArraysOfArrays
using CairoMakie
using IntervalSets
using LegendDSP
using LegendMakie
using RadiationDetectorDSP
using RadiationDetectorSignals
using Random
using StableRNGs
using Statistics
using StatsBase
using Unitful

CairoMakie.activate!(type = "png")
update_theme!(fontsize = 14)

# ## Synthetic Events

# The waveforms have the same properties as in the energy tutorial: a 16 ns sampling period, a decay time of 40 µs,
# a baseline offset, white noise and a random walk of the preamplifier output. The pulse starts at (50 ± 2) µs, so that
# the pre-trigger baseline of at least 48 µs is longer than all filters used below. In contrast to the energy tutorial,
# the width of the current pulse is almost the same for all events (100 ns with a spread of 2%), so that the
# A/E of single-site events has a narrow distribution.
#
# An event has a current pulse that is a Gaussian of width 100 ns. For a multi-site event, a fraction `f` of the energy
# is carried by a second Gaussian pulse that is delayed by `delay`. The fraction is between 15% and 50% and skewed
# towards small values, so that most multi-site events are just below the single-site band. 15% of the multi-site events have
# a third pulse with a fraction of 10% to 40%. Delays between 50 and 480 ns are drawn, so that the pulses can overlap partly. These are
# the delays between the arrival of charges deposited at
# different distances from the contact. A fraction `q` of the charge of an event can also arrive late,
# as an exponential with a time constant of 1.5 µs after the main pulse. This is used in the LQ section.
#
# The energies follow a toy Th-228 calibration spectrum with a continuum that decreases slowly with energy
# and Gaussian lines with the relative intensities and the single-site fractions given in the table.
# The single-site fractions are typical for the full energy peaks and escape peaks of large germanium detectors:
# the double escape peak (DEP) consists mostly of single-site events, while most events in the single escape peak (SEP) and in the full energy
# peak (FEP) are multi-site.
#
# | Population | Energy (keV) | Events | Single-site fraction |
# |:-----------|:-------------|-------:|---------------------:|
# | Compton continuum | 300 - 2600 | 4000 | 0.40 |
# | FEP (Tl-208) | 583.2 | 350 | 0.30 |
# | FEP (Tl-208) | 727.3 | 80 | 0.25 |
# | FEP (Bi-212) | 860.6 | 50 | 0.20 |
# | DEP | 1592.5 | 150 | 0.90 |
# | FEP (Bi-212) | 1620.7 | 35 | 0.12 |
# | SEP | 2103.5 | 80 | 0.05 |
# | FEP | 2614.5 | 500 | 0.08 |

function simulate_aoe_waveforms(rng, energies, fractions, delays, late = zeros(length(energies));
        fractions3 = zeros(length(energies)), delays3 = fill(0.0u"ns", length(energies)), dt = 16.0u"ns", n = 6000,
        t_start = 50.0u"µs", τ = 40.0u"µs", τ_late = 1.5u"µs", gain = 2f0, white_noise = 10f0, walk_noise = 0.07f0)
    time = range(0.0u"µs", step = dt, length = n)
    t = Float32.(ustrip.(u"ns", time))
    step_ns = Float32(ustrip(u"ns", dt))
    τ_late_ns = Float32(ustrip(u"ns", τ_late))
    a = Float32(exp(-ustrip(NoUnits, dt / τ)))
    signal = Matrix{Float32}(undef, n, length(energies))
    pulse = similar(t)
    for (j, E) in enumerate(energies)
        t0 = Float32(ustrip(u"ns", t_start + (rand(rng) - 0.5) * 4.0u"µs"))
        width = Float32(100 * (1 + 0.02 * randn(rng)))
        f = Float32(fractions[j])
        d = Float32(ustrip(u"ns", delays[j]))
        f3 = Float32(fractions3[j])
        d3 = Float32(ustrip(u"ns", delays3[j]))
        q = Float32(late[j])
        @. pulse = ((1 - f - f3) * exp(-0.5f0 * ((t - t0) / width)^2) + f * exp(-0.5f0 * ((t - t0 - d) / width)^2) +
                    f3 * exp(-0.5f0 * ((t - t0 - d3) / width)^2)) * step_ns / (width * sqrt(2f0 * π))
        charge = (1 - q) .* cumsum(pulse)
        @. charge += q * ifelse(t > t0, 1 - exp(-(t - t0) / τ_late_ns), 0f0)
        charge .*= gain * Float32(E)
        offset = 2000f0 + 30f0 * randn(rng, Float32)
        y = 0f0
        previous = 0f0
        for i in eachindex(t)
            y = a * y + (charge[i] - previous) + walk_noise * randn(rng, Float32)
            previous = charge[i]
            signal[i, j] = y + offset + white_noise * randn(rng, Float32)
        end
    end
    return ArrayOfRDWaveforms((fill(time, length(energies)), VectorOfSimilarVectors(signal)))
end

rng = StableRNG(7)
gain = 2.0 # ADC counts per keV
τ = 40.0u"µs"

lines = [583.2, 727.3, 860.6, 1592.5, 1620.7, 2103.5, 2614.5] # keV
n_lines = [350, 80, 50, 150, 35, 80, 500]
p_sse = [0.30, 0.25, 0.20, 0.90, 0.12, 0.05, 0.08]
n_continuum, p_sse_continuum = 4000, 0.40

energies = Float64[]
kind = Int[] # index in `lines`, 0 for the continuum
is_sse = Bool[]
for (i, (E, n, p)) in enumerate(zip(lines, n_lines, p_sse))
    append!(energies, E .+ sqrt(0.2^2 + 1e-4 * E) .* randn(rng, n))
    append!(kind, fill(i, n))
    append!(is_sse, rand(rng, n) .< p)
end
x = rand(rng, 3 * n_continuum)
x = x[rand(rng, length(x)) .< 1 .- 0.6 .* x][1:n_continuum] # density falling linearly to 40% of its start value
append!(energies, 300 .+ 2300 .* x)
append!(kind, fill(0, n_continuum))
append!(is_sse, rand(rng, n_continuum) .< p_sse_continuum)

is_mse = .!is_sse
fractions = [is_mse[i] ? 0.15 + 0.35 * rand(rng)^1.5 : 0.0 for i in eachindex(energies)]
delays = (50 .+ 430 .* rand(rng, length(energies))) .* u"ns"
three_site = is_mse .& (rand(rng, length(energies)) .< 0.15)
fractions3 = [three_site[i] ? 0.1 + 0.3 * rand(rng)^1.5 : 0.0 for i in eachindex(energies)]
delays3 = delays .+ (50 .+ 300 .* rand(rng, length(energies))) .* u"ns"

wvfs = simulate_aoe_waveforms(rng, energies, fractions, delays; fractions3, delays3, gain = Float32(gain), τ = τ)
length(wvfs)

# ## Baseline, Pole-Zero Correction and Energy

# The preprocessing is the same as in the energy tutorial. The energy is determined with a trapezoidal filter of
# 4 µs rise time and 2 µs flat-top time. The pick-off is taken at the center of the flat top, `rise + flat / 2`
# after the time at which the waveform reaches 50% of its height. The energies are converted to keV
# with the gain of the simulation. Real data are calibrated as described in the energy tutorial.

bl = signalstats.(wvfs, 0u"µs", 40u"µs")
wvfs_pz = InvCRFilter(τ).(shift_waveform.(wvfs, -bl.mean));

e_101010 = maximum.(TrapezoidalChargeFilter(10u"µs", 10u"µs").(wvfs_pz).signal)
t50 = get_threshold(wvfs_pz, 0.5 .* e_101010; mintot = 600u"ns")
pickoff = SignalEstimator(PolynomialDNI(3, 800u"ns"))
e = pickoff.(TrapezoidalChargeFilter(4.0u"µs", 2.0u"µs").(wvfs_pz), t50 .+ 5.0u"µs") ./ gain;

# ## Current Pulse

# The current is the time derivative of the charge waveform. The derivative of the noisy waveform is dominated by noise,
# so it has to be smoothed. A Savitzky-Golay filter fits a polynomial of degree 3 in a window of 180 ns around each sample
# and gives the derivative of that polynomial. The argument `1` selects the first derivative.
# The result is in units of waveform samples per sample. The window for the maximum, `current_window`, contains the start of the pulse.

current_window = 40u"µs" .. 70u"µs"
wvfs_sg = SavitzkyGolayFilter(180u"ns", 3, 1).(wvfs_pz);

# The maximum of the current in the window around the start of the pulse is A. [`get_a`](@ref) returns the maximum,
# interpolated with a parabola through the three samples around it, and its time.
# It is based on [`get_wvf_maximum`](@ref).

a_sg = get_a(wvfs_sg, current_window)
keys(a_sg[1])

# The current of a single-site and of a multi-site event with similar energy, with the time axis shifted
# to the time at which the charge reaches 50% of its height:

i_sse = findfirst(i -> kind[i] == 0 && is_sse[i] && 1400 < energies[i] < 1500, eachindex(energies))
i_mse = findfirst(i -> kind[i] == 0 && is_mse[i] && 1400 < energies[i] < 1500 && delays[i] > 200u"ns" && fractions[i] > 0.3 && !three_site[i], eachindex(energies))

fig = Figure(size = (800, 380))
ax = Axis(fig[1, 1], xlabel = "Time relative to t50 (µs)", ylabel = "Current (counts / sample)", limits = (-1.5, 2.5, nothing, nothing))
for (i, label) in ((i_sse, "single-site"), (i_mse, "multi-site"))
    lines!(ax, ustrip.(u"µs", wvfs_sg[i].time) .- ustrip(u"µs", t50[i]), wvfs_sg[i].signal, label = label)
    hlines!(ax, a_sg.max[i], color = :black, linestyle = :dash)
end
axislegend(ax, position = :rt)
fig

# The multi-site event has the lower maximum for a similar energy.

# ## Alternative Current Filters

# LegendDSP provides other filters that can be used to smooth the current:
#
# - [`WeightedSavitzkyGolayFilter`](@ref) fits the polynomial with weights that fall off towards the edge of the window.
#   The weight type 2 is a Hann window. The filter returns the smoothed signal, which is differentiated
#   afterwards with [`DerivativeFilter`](@ref).
# - [`MovingWindowMultiFilter`](@ref) averages the derivative of the waveform in a moving window. Three
#   passes alternating in direction are used, which makes the result symmetric.
#
# The window lengths are of the same order, because all three filters reduce the noise with the same bandwidth.

wvfs_deriv = DerivativeFilter(1).(wvfs_pz)
wvfs_wsg = DerivativeFilter(1).(WeightedSavitzkyGolayFilter(length = 300u"ns", degree = 3, weightType = 2).(wvfs_pz))
wvfs_mwa = MovingWindowMultiFilter(100u"ns"; passes = 3).(wvfs_deriv)

a_wsg = get_a(wvfs_wsg, current_window)
a_mwa = get_a(wvfs_mwa, current_window)
a_raw = get_a(wvfs_deriv, current_window);

# The filters are compared by the width of the A/E distribution of the single-site events in the DEP.
# The normalization is explained in the next section. The DEP events are selected by their energy in a window of
# ±5 keV around the peak. A/E is divided by the median of these events,
# most of which are single-site. The width is the standard deviation of the A/E values (in %)
# estimated robustly from the median absolute deviation of the DEP events above the median,
# which excludes the multi-site tail on the low side.

dep = findall(abs.(e .- 1592.5) .< 5)
upper_width(x) = 100 * mad(x[x .> median(x)] .- median(x), center = 0, normalize = true) / median(x)

widths = map((raw = a_raw, sg = a_sg, wsg = a_wsg, mwa = a_mwa)) do a
    aoe = a.max[dep] ./ e[dep]
    upper_width(aoe)
end

# The width of the unfiltered derivative is dominated by noise. The three smoothing filters give a similar, smaller width.

# ## A/E

# A/E is the ratio of the maximum of the current to the energy. Its absolute scale depends on the gain
# and the units of the filters. A/E is therefore normalized by the median in the DEP, which is dominated by single-site events, so
# that the single-site band is at 1.

aoe = a_sg.max ./ e
aoe ./= median(aoe[dep])

# The histogram of A/E has a logarithmic count axis. Empty bins are set to 0.1, below the lower limit of the axis,
# so that the line drops to the axis for empty bins.

bins = 0.4:0.01:1.2
fig = Figure(size = (800, 380))
ax = Axis(fig[1, 1], xlabel = "A/E (normalized)", ylabel = "Events / 0.01", yscale = log10, limits = (0.4, 1.2, 0.8, nothing))
for (selection, label) in ((eachindex(aoe), "all events"), (dep, "DEP window"))
    counts = fit(Histogram, aoe[selection], bins).weights
    stairs!(ax, bins, [max.(counts, 0.1); 0.1], step = :post, label = label)
end
axislegend(ax, position = :lt)
fig

# The single-site events form a peak at 1. Multi-site events are at lower values
# and form a tail. The tail of the DEP window contains the multi-site events of the peak and of the continuum below it.
# The following figure shows the energy and A/E of all events in a two-dimensional histogram.
# Empty bins are not drawn.

fig = Figure(size = (800, 450))
ax = Axis(fig[1, 1], xlabel = "Energy (keV)", ylabel = "A/E (normalized)", limits = (400, 2700, 0.4, 1.2))
h2 = fit(Histogram, (e, aoe), (400:25:2700, 0.4:0.02:1.2))
hm = heatmap!(ax, h2.edges[1], h2.edges[2], replace(h2.weights, 0 => NaN), colorscale = log10)
Colorbar(fig[1, 2], hm, label = "Events")
fig

# The single-site events form a band at 1 and the multi-site events are spread below it. The DEP at 1592.5 keV
# stands out in the band as a column of mostly single-site events. The single-site fraction of the
# SEP at 2103.5 keV and of the FEP at 2614.5 keV is small, so that both lines are mostly below the band.
# Lines with a large multi-site fraction are visible as columns in the tail.
#
# The calibration of A/E, with a correction of its energy dependence and a cut that is applied
# to select single-site events, is done with LegendSpecFits.jl.
# The fraction of events above a cut, here at 0.9, in windows of ±4 keV around the peaks and in a part of the continuum
# without peaks, is read from the sample directly:

cut = 0.9
survival(lo, hi) = mean(aoe[lo .< e .< hi] .> cut)
(dep = survival(1588.5, 1596.5), sep = survival(2099.5, 2107.5), fep = survival(2610.5, 2618.5), continuum = survival(1700, 2050))

# In this sample, more than 80% of the events in the DEP window survive, and only
# 15% to 20% of the SEP and FEP events. Of the events in the continuum, about half survive.

# ## Charge Drift: Q-Drift

# The charge-drift parameter Q-drift is sensitive to the arrival time of the charge. [`get_qdrift`](@ref)
# integrates the waveform and compares the charge collected in two successive time intervals after t0, the time at which the
# pulse starts: the first from `t0` to `t0 + first(Δt)` and the second from there to `t0 + last(Δt)`.
# With `Δt = 2.5 µs:0.1 µs:5 µs` the two intervals have the same length. For an ideal step starting at t0
# both intervals contain the same charge and the difference is zero. If charge arrives later,
# the second interval has more charge. Q-drift is proportional to the energy and to the delay of the charge relative
# to t0, so Q-drift divided by the energy is a delay.

t0 = get_t0(wvfs_pz, 80; mintot = 600u"ns")
qdrift = get_qdrift(wvfs_pz, t0, 2.5u"µs":0.1u"µs":5.0u"µs");

# The integrator filter used by `get_qdrift` sums samples, so Q-drift has the units of counts times samples.
# Dividing by the energy in counts and multiplying by the sampling period gives a time. For single-site events it is
# the time between t0 and the arrival of the charge, which is a few hundred nanoseconds because t0 is
# the first time the trigger filter crosses its threshold, before the charge reaches 50% of its height.
# For multi-site events the later charge adds to it. The two populations are separated here by the simulated truth:

step_t = step(wvfs[1].time)
qdrift_time = uconvert.(u"ns", qdrift .* step_t ./ (e .* gain))
sel = findall((kind .== 0) .& (e .> 1000))
bins = 100:10:600

fig = Figure(size = (800, 380))
ax = Axis(fig[1, 1], xlabel = "Q-drift time (ns)", ylabel = "Events / 10 ns")
stephist!(ax, ustrip.(u"ns", qdrift_time[sel[is_sse[sel]]]), bins = bins, label = "single-site")
stephist!(ax, ustrip.(u"ns", qdrift_time[sel[is_mse[sel]]]), bins = bins, label = "multi-site")
axislegend(ax, position = :rt)
fig

# ## Late Charge: LQ

# Charge that is collected late, for example charge that drifts slowly close to the surface of the detector, makes the
# waveform rise slowly after the main pulse. LQ is the difference of the areas under the pole-zero corrected waveform in two
# successive windows after the time at which the waveform reaches a fraction of its height. [`get_lq`](@ref)
# takes the integrated waveform, the start time, and the lengths `early`, `gap` and `late` of the early window, of a gap
# and of the late window: the early window is `[t_start, t_start + early]` and the late window starts after the gap.
# The result is the area of the late window minus the area of the early window. For a pulse with a flat top,
# both areas are equal and LQ is close to zero. For late charge, the waveform is still rising
# in the late window and LQ is positive. [`get_qdrift`](@ref) is the special case of the early window
# starting at t0 without a gap.
#
# Here LQ is computed for 600 single-site events with energies between 1000 and 2600 keV, of which 30% have late charge, a fraction
# of 2% to 8% of the energy that arrives with a time constant of 1.5 µs.

n_lq = 600
e_lq_true = 1000 .+ 1600 .* rand(rng, n_lq)
late = [rand(rng) < 0.3 ? 0.02 + 0.06 * rand(rng) : 0.0 for _ in 1:n_lq]
wvfs_lq = simulate_aoe_waveforms(rng, e_lq_true, zeros(n_lq), fill(0.0u"ns", n_lq), late; gain = Float32(gain), τ = τ)

bl_lq = signalstats.(wvfs_lq, 0u"µs", 40u"µs")
wvfs_lq_pz = InvCRFilter(τ).(shift_waveform.(wvfs_lq, -bl_lq.mean))
e_lq = maximum.(TrapezoidalChargeFilter(10u"µs", 10u"µs").(wvfs_lq_pz).signal)
t80_lq = get_threshold(wvfs_lq_pz, 0.8 .* e_lq; mintot = 600u"ns");

# The windows start at t80, the time at which the waveform reaches 80% of its height, because the rise of the
# charge is complete shortly after it. The windows have a length of 2.5 µs each and no gap, which is the configuration
# `lq_windows = [[2.5u"µs", 0.0u"µs", 2.5u"µs"]]` of [`DSPConfig`](@ref) used in the energy tutorial.

early, gap, late_window = 2.5u"µs", 0.0u"µs", 2.5u"µs"
lq = get_lq(IntegratorFilter(1).(wvfs_lq_pz), t80_lq, early, gap, late_window);

# One event without and one with late charge, normalized to their heights. The upper panel shows the areas under the
# waveform in the early window (blue) and in the late window (red), which [`get_lq`](@ref) integrates. LQ is the red area minus the blue
# area. Both areas are close to the height times the window length, so the lower panel shows what makes them different:
# the area between the waveform and the flat-top level, that is the charge that has not arrived yet.
# LQ is the deficit of the early window minus the deficit of the late window. It is larger for the event with late
# charge, because most of its deficit is in the early window and much less is left in the late window.

i_ss = findfirst(i -> late[i] == 0 && 1500 < e_lq_true[i] < 1900, 1:n_lq)
i_lc = findfirst(i -> late[i] > 0.07 && 1500 < e_lq_true[i] < 1900, 1:n_lq)

fig = Figure(size = (800, 640))
ax1 = Axis(fig[1, 1], ylabel = "Signal (fraction of the height)", limits = (-1.5, 7.5, -0.05, 1.15))
ax2 = Axis(fig[2, 1], xlabel = "Time relative to t80 (µs)", ylabel = "Signal (fraction of the height)", limits = (-1.5, 7.5, 0.85, 1.03))
linkxaxes!(ax1, ax2)
for (i, label, color) in ((i_ss, "no late charge", :black), (i_lc, "late charge", :darkorange))
    t = ustrip.(u"µs", wvfs_lq_pz[i].time) .- ustrip(u"µs", t80_lq[i])
    y = wvfs_lq_pz[i].signal ./ e_lq[i]
    in_early = 0 .<= t .<= ustrip(u"µs", early)
    in_late = ustrip(u"µs", early + gap) .<= t .<= ustrip(u"µs", early + gap + late_window)
    for ax in (ax1, ax2)
        lines!(ax, t, y, color = color, label = label)
    end
    i == i_lc || continue
    band!(ax1, t[in_early], zeros(count(in_early)), y[in_early], color = (:royalblue, 0.25))
    band!(ax1, t[in_late], zeros(count(in_late)), y[in_late], color = (:crimson, 0.25))
    band!(ax2, t[in_early], y[in_early], ones(count(in_early)), color = (:royalblue, 0.4))
    band!(ax2, t[in_late], y[in_late], ones(count(in_late)), color = (:crimson, 0.4))
end
hlines!(ax2, 1.0, color = :gray, linestyle = :dash)
axislegend(ax2, position = :rb)
fig

# LQ is proportional to the energy. Dividing it by the energy in counts and multiplying with the sampling period
# gives a time, as for Q-drift. The two populations are separated in LQ over the energy range:

lq_time = uconvert.(u"ns", lq .* step_t ./ e_lq)
has_late = late .> 0

fig = Figure(size = (800, 380))
ax = Axis(fig[1, 1], xlabel = "Energy (keV)", ylabel = "LQ / E (ns)", limits = (nothing, nothing, 0, 110))
scatter!(ax, e_lq[.!has_late] ./ gain, ustrip.(u"ns", lq_time[.!has_late]), markersize = 5, label = "no late charge")
scatter!(ax, e_lq[has_late] ./ gain, ustrip.(u"ns", lq_time[has_late]), markersize = 5, label = "late charge")
axislegend(ax, position = :lt)
fig

# Events without late charge scatter around a small value that is given by the shape of the rise. The late-charge events
# are spread to larger values, according to the simulated fraction of late charge. In an analysis, a cut is set on the
# energy-normalized LQ, as for A/E. [`dsp_icpc`](@ref) returns LQ for every window of `config.lq_windows` and for the start times
# `t80`, `t90` and `t_amax`, in the columns `lq80_<early><gap><late>`, `lq90_<early><gap><late>` and `lqamax_<early><gap><late>`.
# The table of the energy tutorial contains, for example, the column `lq80_250025`.
