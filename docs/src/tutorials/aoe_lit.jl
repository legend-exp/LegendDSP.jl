# # A/E Pulse-Shape Parameter

# Single-site events (SSE) deposit their energy in a small volume of the detector. Their current pulse is short and
# its maximum A is proportional to the energy E. Multi-site events (MSE), in which the energy is deposited
# at several places, produce a sum of pulses that arrive at different times. The maximum of the summed
# current is lower for the same energy. The ratio A/E therefore separates SSE from MSE.
# This tutorial computes A/E from simulated waveforms and compares alternative current filters.
# The charge-drift parameter Q-drift is computed as well.

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
# a baseline offset, white noise and a random walk of the preamplifier output. In contrast to that tutorial,
# the width of the current pulse is almost the same for all events (100 ns with a spread of 2%), so that the
# A/E of single-site events has a narrow distribution.
#
# An event has a current pulse that is a Gaussian of width 100 ns. For a multi-site event, a fraction `f` of the energy
# is carried by a second Gaussian pulse that is delayed by `delay`. Delays of tens
# to hundreds of nanoseconds are drawn, which are the delays between the arrival of charges deposited at
# different distances from the contact.
#
# Two populations are simulated:
#
# - Events of the double escape peak (DEP) at 1592.5 keV (Gaussian with a standard deviation of 0.6 keV),
#   of which 10% are multi-site.
# - Compton continuum events with a flat spectrum between 400 and 2000 keV, of which 70% are multi-site.

function simulate_aoe_waveforms(rng, energies, fractions, delays; dt = 16.0u"ns", n = 5000, t_start = 25.0u"µs",
        τ = 40.0u"µs", gain = 2f0, white_noise = 10f0, walk_noise = 0.07f0)
    time = range(0.0u"µs", step = dt, length = n)
    t = Float32.(ustrip.(u"ns", time))
    step_ns = Float32(ustrip(u"ns", dt))
    a = Float32(exp(-ustrip(NoUnits, dt / τ)))
    signal = Matrix{Float32}(undef, n, length(energies))
    pulse = similar(t)
    for (j, E) in enumerate(energies)
        t0 = Float32(ustrip(u"ns", t_start + (rand(rng) - 0.5) * 4.0u"µs"))
        width = Float32(100 * (1 + 0.02 * randn(rng)))
        f = Float32(fractions[j])
        d = Float32(ustrip(u"ns", delays[j]))
        @. pulse = ((1 - f) * exp(-0.5f0 * ((t - t0) / width)^2) + f * exp(-0.5f0 * ((t - t0 - d) / width)^2)) *
                   step_ns / (width * sqrt(2f0 * π))
        charge = cumsum(pulse)
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

n_dep, n_compton = 600, 6000
energies = [1592.5 .+ 0.6 .* randn(rng, n_dep); 400 .+ 1600 .* rand(rng, n_compton)]
is_dep = [trues(n_dep); falses(n_compton)]
is_mse = [rand(rng, n_dep) .< 0.1; rand(rng, n_compton) .< 0.7]
fractions = [is_mse[i] ? 0.2 + 0.6 * rand(rng) : 0.0 for i in eachindex(energies)]
delays = (50 .+ 350 .* rand(rng, length(energies))) .* u"ns"

wvfs = simulate_aoe_waveforms(rng, energies, fractions, delays; gain = Float32(gain), τ = τ)
length(wvfs)

# ## Baseline, Pole-Zero Correction and Energy

# The preprocessing is the same as in the energy tutorial. The energy is determined with a trapezoidal filter of
# 4 µs rise time and 2 µs flat-top time. The pick-off is taken at the center of the flat top, `rise + flat / 2`
# after the time at which the waveform reaches 50% of its height. The energies are converted to keV
# with the gain of the simulation. Real data are calibrated as described in the energy tutorial.

bl = signalstats.(wvfs, 0u"µs", 20u"µs")
wvfs_pz = InvCRFilter(τ).(shift_waveform.(wvfs, -bl.mean));

e_101010 = maximum.(TrapezoidalChargeFilter(10u"µs", 10u"µs").(wvfs_pz).signal)
t50 = get_threshold(wvfs_pz, 0.5 .* e_101010; mintot = 600u"ns")
pickoff = SignalEstimator(PolynomialDNI(3, 800u"ns"))
e = pickoff.(TrapezoidalChargeFilter(4.0u"µs", 2.0u"µs").(wvfs_pz), t50 .+ 5.0u"µs") ./ gain;

# ## Current Pulse

# The current is the time derivative of the charge waveform. The derivative of the noisy waveform is dominated by noise,
# so it has to be smoothed. A Savitzky-Golay filter fits a polynomial of degree 3 in a window of 180 ns around each sample
# and gives the derivative of that polynomial. The argument `1` selects the first derivative.
# The result is in units of waveform samples per sample.

current_window = 15u"µs" .. 45u"µs"
wvfs_sg = SavitzkyGolayFilter(180u"ns", 3, 1).(wvfs_pz);

# The maximum of the current in the window around the start of the pulse is A. [`get_a`](@ref) returns the maximum,
# interpolated with a parabola through the three samples around it, and its time.
# It is based on [`get_wvf_maximum`](@ref).

a_sg = get_a(wvfs_sg, current_window)
keys(a_sg[1])

# The current of a single-site and of a multi-site event with similar energy, with the time axis shifted
# to the time at which the charge reaches 50% of its height:

i_sse = findfirst(i -> !is_mse[i] && !is_dep[i] && 1400 < energies[i] < 1500, eachindex(energies))
i_mse = findfirst(i -> is_mse[i] && !is_dep[i] && 1400 < energies[i] < 1500 && delays[i] > 300u"ns" && fractions[i] > 0.4, eachindex(energies))

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
# The normalization is explained in the next section. A/E is divided by the median of the DEP events,
# most of which are single-site. The width is the standard deviation of the A/E values (in %)
# estimated robustly from the median absolute deviation of the DEP events above the median,
# which excludes the multi-site tail on the low side.

dep = findall(is_dep)
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

fig = Figure(size = (800, 380))
ax = Axis(fig[1, 1], xlabel = "A/E (normalized)", ylabel = "Events / 0.01", yscale = log10, limits = (0.4, 1.2, 0.5, nothing))
stephist!(ax, aoe, bins = 0.4:0.01:1.2, label = "all events")
stephist!(ax, aoe[dep], bins = 0.4:0.01:1.2, label = "DEP")
axislegend(ax, position = :lt)
fig

# The single-site events form a peak at 1. Multi-site events are at lower values
# and form a tail. The tail of the DEP is the 10% of multi-site events in it.
# The following figure shows the energy and A/E of all events in a two-dimensional histogram.

fig = Figure(size = (800, 450))
ax = Axis(fig[1, 1], xlabel = "Energy (keV)", ylabel = "A/E (normalized)", limits = (400, 2000, 0.4, 1.2))
h2 = fit(Histogram, (e, aoe), (400:40:2000, 0.4:0.02:1.2))
counts = replace(h2.weights, 0 => NaN)
hm = heatmap!(ax, h2.edges[1], h2.edges[2], counts, colorscale = log10)
Colorbar(fig[1, 2], hm, label = "Events")
fig

# The single-site events form a band at 1 and the multi-site events are spread below it. The DEP at 1592.5 keV is
# visible as a column of mostly single-site events.
#
# The calibration of A/E, with a correction of its energy dependence and a cut that is applied
# to select single-site events, is done with LegendSpecFits.jl.
# For the DEP events the fraction below a cut can be read from the sample directly:

cut = 0.9
(survival_dep = mean(aoe[dep] .> cut), survival_continuum = mean(aoe[.!is_dep] .> cut))

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
sel = findall(.!is_dep .& (e .> 1000))
bins = 100:10:450

fig = Figure(size = (800, 380))
ax = Axis(fig[1, 1], xlabel = "Q-drift time (ns)", ylabel = "Events / 10 ns")
stephist!(ax, ustrip.(u"ns", qdrift_time[sel[.!is_mse[sel]]]), bins = bins, label = "single-site")
stephist!(ax, ustrip.(u"ns", qdrift_time[sel[is_mse[sel]]]), bins = bins, label = "multi-site")
axislegend(ax, position = :rt)
fig

# The function [`get_lq`](@ref) and the filter `IntegratorFilter` from RadiationDetectorDSP give a related late-charge
# parameter, see [Waveform Parameters](@ref). [`dsp_icpc`](@ref) returns A/E ingredients (`a_sg`, `a_mwa`, `a_raw`),
# `qdrift` and the LQ parameters for every event, as shown at the end of the [Energy Reconstruction](@ref) tutorial.
