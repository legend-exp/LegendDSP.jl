# # Energy Reconstruction

# This tutorial reconstructs the energy of charge waveforms from a high-purity germanium (HPGe) detector,
# step by step: baseline subtraction, decay-time estimation, pole-zero correction, determination of the
# signal start, energy filters and a final energy spectrum. At the end, the same steps are carried out in one call
# with [`dsp_icpc`](@ref).

# The waveforms are simulated, so that the example runs without detector data and the
# results of each step can be compared with the values used in the simulation.
# For the reading of real data see the documentation of
# [LegendDataManagement.jl](https://legend-exp.github.io/LegendDataManagement.jl/stable/).
# The manual pages [HPGe DSP Chains](@ref), [Filters](@ref), [Waveform Parameters](@ref) and
# [Decay Time](@ref) describe the functions used here in more detail.

# ## Packages

using CairoMakie
using IntervalSets
using LegendDSP
using LegendMakie
using PropDicts
using RadiationDetectorDSP
using RadiationDetectorSignals
using Random
using StableRNGs
using Statistics
using StatsBase
using TypedTables
using Unitful
using ArraysOfArrays

CairoMakie.activate!(type = "png")
update_theme!(fontsize = 14)

# ## Synthetic Waveforms

# A charge-sensitive preamplifier turns the charge collected in the detector into a step that
# decays exponentially with a time constant `τ`. The simulated waveforms have the following properties:
#
# - The sampling period is 16 ns and a trace has 5000 samples (80 µs).
# - The charge rises with the shape of the integral of a Gaussian current pulse, with a width between 80 and 140 ns
#   for each event. The 10%-90% rise time is therefore between 200 and 360 ns.
#   The pulse starts at (25 ± 2) µs.
# - The decay time is `τ = 40 µs`.
# - Every waveform has an individual baseline offset of 2000 ± 30 ADC counts. White noise (standard deviation of 10 counts)
#   is added to every sample. A random walk of the preamplifier output, which models the parallel noise of the
#   input stage, is added before the decay.
# - The gain is 2 ADC counts per keV.
#
# The energies follow a toy Th-228 calibration spectrum. It contains Gaussian lines at 583.2 keV, 860.6 keV, the
# double escape peak (DEP) at 1592.5 keV, the single escape peak (SEP) at 2103.5 keV and the full energy peak (FEP) at 2614.5 keV,
# and a flat continuum between 150 and 2614 keV.
# The intrinsic standard deviation of a line is `sqrt(0.2² + 10⁻⁴ E)` in keV for the energy `E` in keV.

function simulate_charge_waveforms(rng, energies; dt = 16.0u"ns", n = 5000, t_start = 25.0u"µs", τ = 40.0u"µs",
        gain = 2f0, white_noise = 10f0, walk_noise = 0.07f0)
    time = range(0.0u"µs", step = dt, length = n)
    t = Float32.(ustrip.(u"ns", time))
    a = Float32(exp(-ustrip(NoUnits, dt / τ)))
    signal = Matrix{Float32}(undef, n, length(energies))
    t0 = similar(energies, typeof(t_start))
    baseline = similar(energies, Float32)
    pulse = similar(t)
    for (j, E) in enumerate(energies)
        t0[j] = t_start + (rand(rng) - 0.5) * 4.0u"µs"
        width = Float32(80 + 60 * rand(rng))
        @. pulse = exp(-0.5f0 * ((t - Float32(ustrip(u"ns", t0[j]))) / width)^2)
        charge = cumsum(pulse)
        charge .*= gain * Float32(E) * Float32(ustrip(u"ns", dt)) / (width * sqrt(2f0 * π))
        baseline[j] = 2000f0 + 30f0 * randn(rng, Float32)
        y = 0f0
        previous = 0f0
        for i in eachindex(t)
            y = a * y + (charge[i] - previous) + walk_noise * randn(rng, Float32)
            previous = charge[i]
            signal[i, j] = y + baseline[j] + white_noise * randn(rng, Float32)
        end
    end
    wvfs = ArrayOfRDWaveforms((fill(time, length(energies)), VectorOfSimilarVectors(signal)))
    return (; wvfs, t0, baseline, energies)
end

# The container `ArrayOfRDWaveforms` stores the waveforms of many events in two arrays, one for the time axes
# and one for the signals. Indexing it returns single `RDWaveform`s. All
# functions of LegendDSP and RadiationDetectorDSP can be broadcast over it.

rng = StableRNG(42)
gain = 2.0 # ADC counts per keV
τ_true = 40.0u"µs"

lines = [583.2, 860.6, 1592.5, 2103.5, 2614.5] # keV
n_lines = [400, 150, 250, 250, 800]

energies = Float64[]
line_index = Int[] # 0 for continuum events
for (i, (E, n)) in enumerate(zip(lines, n_lines))
    append!(energies, E .+ sqrt(0.2^2 + 1e-4 * E) .* randn(rng, n))
    append!(line_index, fill(i, n))
end
append!(energies, 150 .+ (2614 - 150) .* rand(rng, 2000))
append!(line_index, fill(0, 2000))
order = randperm(rng, length(energies))
energies, line_index = energies[order], line_index[order]

sim = simulate_charge_waveforms(rng, energies; gain = Float32(gain), τ = τ_true)
wvfs = sim.wvfs
length(wvfs)

# Three waveforms with different energies:

fig = Figure(size = (800, 350))
ax = Axis(fig[1, 1], xlabel = "Time (µs)", ylabel = "Signal (ADC counts)")
for i in (findfirst(==(5), line_index), findfirst(==(1), line_index), findfirst(==(0), line_index))
    LegendMakie.waveformplot!(ax, wvfs[i], label = "$(round(energies[i], digits = 1)) keV")
end
axislegend(ax, position = :rt)
fig

# ## Baseline

# The baseline of each waveform is the mean of the samples before the signal starts.
# `signalstats` from RadiationDetectorDSP returns the mean, standard deviation,
# slope and offset of a linear fit in a time interval. Here the interval from 0 to 20 µs is used.
# `shift_waveform` subtracts the baseline.

bl_window = 0u"µs" .. 20u"µs"
bl = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))
wvfs_bl = shift_waveform.(wvfs, -bl.mean);

# The estimate is compared with the baseline offsets of the simulation. The typical deviation is about one ADC count,
# which is caused by the white noise and mainly by the random walk within the 20 µs window. The example stops
# if the deviation is much larger:

Δbl = bl.mean .- sim.baseline
median(abs.(Δbl)) < 3 || error("baseline estimate deviates from the simulated baseline")
median(abs.(Δbl))

# ## Estimating the Decay Time

# The tail of the waveform decays as `exp(-t/τ)`. `tailstats` fits a straight line to the logarithm of the signal in a
# time interval and returns the decay time `τ`. It needs a baseline-subtracted waveform with positive samples in the
# interval. The tail window is placed after the pulse, from 45 to 78 µs. The logarithm of low-amplitude
# waveforms is dominated by noise, so only events above 1500 keV are used (selected here with the maximum of the waveform).

tail_window = 45u"µs" .. 78u"µs"
tail = tailstats.(wvfs_bl, leftendpoint(tail_window), rightendpoint(tail_window))
high = findall(maximum.(wvfs_bl.signal) .> 1500 * gain)
τ_events = uconvert.(u"µs", tail.τ[high])

fig = Figure(size = (600, 330))
ax = Axis(fig[1, 1], xlabel = "Decay time (µs)", ylabel = "Events")
hist!(ax, ustrip.(u"µs", τ_events), bins = 30)
vlines!(ax, ustrip(u"µs", τ_true), color = :red, label = "simulated value")
axislegend(ax, position = :rt)
fig

# The median of the distribution is used as the decay time of the detector. It agrees with the simulated value:

τ = median(τ_events)
(estimate = τ, simulated = τ_true)

# [`dsp_decay_times`](@ref) does the same for a set of waveforms in one call:

τ_events_dsp = dsp_decay_times(wvfs[high], bl_window, tail_window)
τ_events_dsp ≈ τ_events

# ## Pole-Zero Correction

# The pole-zero correction removes the exponential decay, so that the pulse ends in a flat top at the
# height of the step. `InvCRFilter` inverts the response of the preamplifier for the decay time `τ`.

wvfs_pz = InvCRFilter(τ).(wvfs_bl);

i = findfirst(==(5), line_index)
fig = Figure(size = (800, 350))
ax = Axis(fig[1, 1], xlabel = "Time (µs)", ylabel = "Signal (ADC counts)")
LegendMakie.waveformplot!(ax, wvfs_bl[i], label = "baseline subtracted")
LegendMakie.waveformplot!(ax, wvfs_pz[i], label = "pole-zero corrected")
axislegend(ax, position = :rt)
fig

# The slope of a linear fit in the tail window, in ADC counts per µs, measures the quality of the correction:

slope_before = signalstats.(wvfs_bl, leftendpoint(tail_window), rightendpoint(tail_window)).slope
slope_after = signalstats.(wvfs_pz, leftendpoint(tail_window), rightendpoint(tail_window)).slope
(before = median(slope_before[high]), after = median(slope_after[high]))

# ## Signal Start and Rise Time

# The time of the start of the pulse, t0, is the first time at which a fast asymmetric trapezoidal filter crosses a threshold.
# The filter averages the signal for 40 ns after the point of interest and for 2000 ns before it, with a gap of
# 100 ns, and is therefore a sensitive detector of a steep rise. [`get_t0`](@ref) applies the filter and
# `Intersect` from RadiationDetectorDSP with a minimum time over threshold of `mintot`.
# The threshold is given in ADC counts.

t0 = get_t0(wvfs_pz, 80; mintot = 600u"ns");

# The time of the waveform at 50% of its height, t50, requires a measure of the height that is not affected by
# noise on the flat top. The maximum of a trapezoidal filter with 10 µs rise time and 10 µs flat top is used.
# [`get_threshold`](@ref) intersects the waveform with the threshold of every event.

flt_101010 = TrapezoidalChargeFilter(10u"µs", 10u"µs")
e_101010 = maximum.(flt_101010.(wvfs_pz).signal)
t50 = get_threshold(wvfs_pz, 0.5 .* e_101010; mintot = 600u"ns");

# The simulated pulse is centered at the true start time, where the charge is at 50% of its final value:

(t50_bias = median(t50 .- sim.t0), t50_spread = mad(ustrip.(u"µs", t50 .- sim.t0), normalize = true) * u"µs")

#

i = findfirst(==(5), line_index)
fig = Figure(size = (800, 350))
ax = Axis(fig[1, 1], xlabel = "Time (µs)", ylabel = "Signal (ADC counts)", limits = (22, 29, nothing, nothing))
LegendMakie.waveformplot!(ax, wvfs_pz[i], label = "pole-zero corrected")
vlines!(ax, ustrip(u"µs", t0[i]), color = :green, label = "t0")
vlines!(ax, ustrip(u"µs", t50[i]), color = :red, label = "t50")
hlines!(ax, 0.5 * e_101010[i], color = :red, linestyle = :dash)
axislegend(ax, position = :rb)
fig

# t0 is earlier than t50 because the threshold of the trigger filter is a small fraction of the pulse height.

# ## Energy Filters

# The energy is the height of the step. Filtering the waveform with a finite impulse response filter reduces the noise
# and gives a pulse with a flat top at a height that is proportional to the energy. The flat top makes the result
# insensitive to the exact position of the pick-off.
#
# - `TrapezoidalChargeFilter(rise, flat)`: the difference between two moving averages of the length
#   `rise` separated by `flat`. A step is transformed into a trapezoid.
# - `CUSPChargeFilter` and `ZACChargeFilter`: filters with a cusp-like shape, which are the optimal filters
#   for white series noise. The arguments are the shaping time, the length of the flat top, a decay time of the
#   filter, the total length of the filter and a scale factor. The decay time is set to a very large value, because
#   the pole-zero correction has already been applied. The scale factor is the length of the filter in samples.
#
# The energy is the value of the filtered waveform at a fixed time after t50. The pick-off is
# done with a `SignalEstimator`, which fits a polynomial of degree 3 in a window of 800 ns around the pick-off time.
# For the trapezoid, the pick-off is at the center of the flat top: `t50 + rise + flat / 2`. For the CUSP and ZAC filters,
# it is at half of the filter length.

filter_length = 30.0u"µs"
scale = ustrip(NoUnits, filter_length / step(wvfs[1].time))
flt_trap = TrapezoidalChargeFilter(4.0u"µs", 2.0u"µs")
flt_cusp = CUSPChargeFilter(2.0u"µs", 2.0u"µs", 1e7u"µs", filter_length, scale)
flt_zac = ZACChargeFilter(4.0u"µs", 2.0u"µs", 1e7u"µs", filter_length, scale)
pickoff = SignalEstimator(PolynomialDNI(3, 800u"ns"))

wvfs_trap = flt_trap.(wvfs_pz)
wvfs_cusp = flt_cusp.(wvfs_pz)
wvfs_zac = flt_zac.(wvfs_pz)

e_trap = pickoff.(wvfs_trap, t50 .+ (4u"µs" + 2u"µs" / 2))
e_cusp = pickoff.(wvfs_cusp, t50 .+ filter_length / 2)
e_zac = pickoff.(wvfs_zac, t50 .+ filter_length / 2);

# The three filters applied to the same waveform, with the pick-off times. Note that a filter shortens the time axis.

i = findfirst(==(5), line_index)
fig = Figure(size = (800, 380))
ax = Axis(fig[1, 1], xlabel = "Time (µs)", ylabel = "Signal (ADC counts)")
LegendMakie.waveformplot!(ax, wvfs_pz[i], label = "pole-zero corrected", color = :gray)
LegendMakie.waveformplot!(ax, wvfs_trap[i], label = "trapezoid")
LegendMakie.waveformplot!(ax, wvfs_cusp[i], label = "CUSP")
LegendMakie.waveformplot!(ax, wvfs_zac[i], label = "ZAC")
vlines!(ax, ustrip(u"µs", t50[i] + 5u"µs"), color = :black, linestyle = :dash)
vlines!(ax, ustrip(u"µs", t50[i] + filter_length / 2), color = :black, linestyle = :dot)
axislegend(ax, position = :lt)
fig

# ## Filter Parameters

# The rise time determines the noise of the energy estimate. A long rise time averages more of the white noise,
# but includes more of the random walk. The resolution of the full energy peak is therefore best for an intermediate rise
# time. The full width at half maximum (FWHM) of the peak is estimated from the median absolute deviation of the energies of the events
# in the line, with the nominal gain. The simulation knows which events belong to the full energy peak.

fwhm(e) = 2 * sqrt(2 * log(2)) * mad(e, normalize = true)
fep = findall(==(5), line_index)

flat = 2.0u"µs"
rise_times = [1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 14.0, 20.0] .* u"µs"
fwhm_trap = map(rise_times) do rise
    wvfs_flt = TrapezoidalChargeFilter(rise, flat).(wvfs_pz[fep])
    fwhm(pickoff.(wvfs_flt, t50[fep] .+ (rise + flat / 2)) ./ gain)
end

fig = Figure(size = (600, 330))
ax = Axis(fig[1, 1], xlabel = "Rise time (µs)", ylabel = "FWHM of the FEP (keV)")
scatterlines!(ax, ustrip.(u"µs", rise_times), fwhm_trap)
fig

# The rise time with the best resolution is a few microseconds. In the
# simulation the white noise dominates at short and the random walk at long rise times.
# In production, the optimal filter parameters are found with the grid scans described in [Filter Optimization](@ref).
#
# Resolution of the three filters at the parameters chosen above:

(trap = fwhm(e_trap[fep] ./ gain), cusp = fwhm(e_cusp[fep] ./ gain), zac = fwhm(e_zac[fep] ./ gain))

# ## Energy Spectrum

# The uncalibrated energies are histogrammed. A linear calibration through the peaks is obtained
# by locating each line with the median of the events in a window around the position expected from the nominal gain
# and a least-squares fit of the positions as a function of the literature energies.
# The simulated gain is exactly known here; the fit recovers it up to small deviations. The windows contain a small continuum, which
# is symmetric to first order.

centroids = map(lines) do E
    window = (E - 8) * gain .. (E + 8) * gain
    median(filter(in(window), e_trap))
end

A = [ones(length(lines)) lines]
offset, slope = A \ centroids
e_cal = (e_trap .- offset) ./ slope

(offset = offset, gain = slope, residuals_keV = (centroids .- offset) ./ slope .- lines)

#

fig = Figure(size = (800, 380))
ax = Axis(fig[1, 1], xlabel = "Energy (keV)", ylabel = "Events / 4 keV", yscale = log10, limits = (0, 2800, 0.5, nothing))
stephist!(ax, e_cal, bins = 0:4:2800)
for E in lines
    vlines!(ax, E, color = :red, linestyle = :dot)
end
fig

# ## One Call with `dsp_icpc`

# [`dsp_icpc`](@ref) applies all of the steps above and many more to a table of waveforms and returns a table with one
# row per event. It takes a [`DSPConfig`](@ref), the decay time `τ`, and a `PropDict` with the optimized filter parameters.
# In production, the configuration comes from the metadata with `DSPConfig(dataprod_config(l200).dsp(filekey).default)`.
# Here it is assembled by hand with the values used in this tutorial. The input table needs
# the columns `waveform`, `baseline`, `timestamp`, `eventnumber` and `daqenergy`, of which only the first is used
# by the steps shown here.

grid(start, stop, step) = Dict(:start => start, :stop => stop, :step => step)
rise_flat(rt, ft) = Dict(:rt => rt * u"µs", :ft => ft * u"µs")

config = DSPConfig(PropDict(Dict(
    :enc_pickoff_trap => 12.0u"µs", :enc_pickoff_zac => 15.0u"µs", :enc_pickoff_cusp => 15.0u"µs",
    :flt_length_cusp => filter_length, :flt_length_zac => filter_length,
    :t0_threshold => 80.0, :inTraceCut_std_threshold => 5.0, :sg_flt_degree => 3,
    :bl_window => Dict(:min => 0.0u"µs", :max => 20.0u"µs"),
    :tail_window => Dict(:min => 45.0u"µs", :max => 78.0u"µs"),
    :current_window => Dict(:min => 15.0u"µs", :max => 45.0u"µs"),
    :auxpz1_window => Dict(:min => 45.0u"µs", :max => 60.0u"µs"),
    :auxpz2_window => Dict(:min => 60.0u"µs", :max => 78.0u"µs"),
    :qdrift_int_length => [2.5u"µs", 5.0u"µs"],
    :lq_windows => [[2.5u"µs", 0.0u"µs", 2.5u"µs"]],
    :e_grid_trap => Dict(:rt => grid(1.0u"µs", 10.0u"µs", 1.0u"µs"), :ft => grid(1.0u"µs", 4.0u"µs", 1.0u"µs")),
    :e_grid_zac => Dict(:rt => grid(1.0u"µs", 10.0u"µs", 1.0u"µs"), :ft => grid(1.0u"µs", 4.0u"µs", 1.0u"µs")),
    :e_grid_cusp => Dict(:rt => grid(1.0u"µs", 10.0u"µs", 1.0u"µs"), :ft => grid(1.0u"µs", 4.0u"µs", 1.0u"µs")),
    :a_grid_wl_sg => grid(100.0u"ns", 500.0u"ns", 50.0u"ns"),
    :a_mwa_windows => [100.0u"ns"], :a_mwa_passes => 3, :a_sg_windows => [200.0u"ns"],
    :flt_defaults => Dict(:sg => 200.0u"ns", :mwa => 100.0u"ns",
        :trap => rise_flat(4.0, 2.0), :zac => rise_flat(4.0, 2.0), :cusp => rise_flat(2.0, 2.0)),
    :kwargs_pars => Dict(
        :sig_interpolation_length => 800.0u"ns", :sig_interpolation_order => 3,
        :int_interpolation_length => 100.0u"ns", :int_interpolation_order => 3,
        :fc_bit_depth => 16, :t0_flt_pars => [40.0u"ns", 100.0u"ns", 2000.0u"ns"],
        :t0_mintot => 600.0u"ns", :tx_mintot => 600.0u"ns", :intrace_mintot => 2.0u"µs"),
)));

# The trapezoid is used with 4 µs rise and 2 µs flat-top time, the CUSP and ZAC filters with the defaults of
# the configuration. An empty `PropDict` as filter parameters selects the defaults in `flt_defaults`.
# The decay time that was estimated above is passed on.

data = Table(
    waveform = wvfs,
    baseline = zeros(Int, length(wvfs)),
    timestamp = fill(0.0u"s", length(wvfs)),
    eventnumber = collect(eachindex(wvfs)),
    daqenergy = zeros(length(wvfs)),
)
dsp = dsp_icpc(data, config, τ, PropDict());

# The table has one column for each parameter. The first ten are:

collect(columnnames(dsp))[1:10]

# The trapezoidal energy `e_trap` is the same as in the step-by-step analysis, because the table uses the same decay
# time, t50 definition and pick-off. The largest difference in keV is:

maximum(abs.(dsp.e_trap .- e_trap)) / gain

# The A/E parameters of the table are the topic of the next tutorial, [A/E Pulse-Shape Parameter](@ref).
