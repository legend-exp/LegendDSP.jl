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

fig = Figure(size = (800, 350))
ax = Axis(fig[1, 1], xlabel = "Time (µs)", ylabel = "Signal (ADC counts)")
for i in (findfirst(==(5), line_index), findfirst(==(1), line_index), findfirst(==(0), line_index))
    LegendMakie.waveformplot!(ax, wvfs[i], label = "$(round(energies[i], digits = 1)) keV")
end
axislegend(ax, position = :rt)
fig

bl_window = 0u"µs" .. 20u"µs"
bl = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))
wvfs_bl = shift_waveform.(wvfs, -bl.mean);

Δbl = bl.mean .- sim.baseline
median(abs.(Δbl)) < 3 || error("baseline estimate deviates from the simulated baseline")
median(abs.(Δbl))

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

τ = median(τ_events)
(estimate = τ, simulated = τ_true)

τ_events_dsp = dsp_decay_times(wvfs[high], bl_window, tail_window)
τ_events_dsp ≈ τ_events

wvfs_pz = InvCRFilter(τ).(wvfs_bl);

i = findfirst(==(5), line_index)
fig = Figure(size = (800, 350))
ax = Axis(fig[1, 1], xlabel = "Time (µs)", ylabel = "Signal (ADC counts)")
LegendMakie.waveformplot!(ax, wvfs_bl[i], label = "baseline subtracted")
LegendMakie.waveformplot!(ax, wvfs_pz[i], label = "pole-zero corrected")
axislegend(ax, position = :rt)
fig

slope_before = signalstats.(wvfs_bl, leftendpoint(tail_window), rightendpoint(tail_window)).slope
slope_after = signalstats.(wvfs_pz, leftendpoint(tail_window), rightendpoint(tail_window)).slope
(before = median(slope_before[high]), after = median(slope_after[high]))

t0 = get_t0(wvfs_pz, 80; mintot = 600u"ns");

flt_101010 = TrapezoidalChargeFilter(10u"µs", 10u"µs")
e_101010 = maximum.(flt_101010.(wvfs_pz).signal)
t50 = get_threshold(wvfs_pz, 0.5 .* e_101010; mintot = 600u"ns");

(t50_bias = median(t50 .- sim.t0), t50_spread = mad(ustrip.(u"µs", t50 .- sim.t0), normalize = true) * u"µs")

i = findfirst(==(5), line_index)
fig = Figure(size = (800, 350))
ax = Axis(fig[1, 1], xlabel = "Time (µs)", ylabel = "Signal (ADC counts)", limits = (22, 29, nothing, nothing))
LegendMakie.waveformplot!(ax, wvfs_pz[i], label = "pole-zero corrected")
vlines!(ax, ustrip(u"µs", t0[i]), color = :green, label = "t0")
vlines!(ax, ustrip(u"µs", t50[i]), color = :red, label = "t50")
hlines!(ax, 0.5 * e_101010[i], color = :red, linestyle = :dash)
axislegend(ax, position = :rb)
fig

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

(trap = fwhm(e_trap[fep] ./ gain), cusp = fwhm(e_cusp[fep] ./ gain), zac = fwhm(e_zac[fep] ./ gain))

centroids = map(lines) do E
    window = (E - 8) * gain .. (E + 8) * gain
    median(filter(in(window), e_trap))
end

A = [ones(length(lines)) lines]
offset, slope = A \ centroids
e_cal = (e_trap .- offset) ./ slope

(offset = offset, gain = slope, residuals_keV = (centroids .- offset) ./ slope .- lines)

fig = Figure(size = (800, 380))
ax = Axis(fig[1, 1], xlabel = "Energy (keV)", ylabel = "Events / 4 keV", yscale = log10, limits = (0, 2800, 0.5, nothing))
stephist!(ax, e_cal, bins = 0:4:2800)
for E in lines
    vlines!(ax, E, color = :red, linestyle = :dot)
end
fig

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

data = Table(
    waveform = wvfs,
    baseline = zeros(Int, length(wvfs)),
    timestamp = fill(0.0u"s", length(wvfs)),
    eventnumber = collect(eachindex(wvfs)),
    daqenergy = zeros(length(wvfs)),
)
dsp = dsp_icpc(data, config, τ, PropDict());

collect(columnnames(dsp))[1:10]

maximum(abs.(dsp.e_trap .- e_trap)) / gain
