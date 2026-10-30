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

bl = signalstats.(wvfs, 0u"µs", 20u"µs")
wvfs_pz = InvCRFilter(τ).(shift_waveform.(wvfs, -bl.mean));

e_101010 = maximum.(TrapezoidalChargeFilter(10u"µs", 10u"µs").(wvfs_pz).signal)
t50 = get_threshold(wvfs_pz, 0.5 .* e_101010; mintot = 600u"ns")
pickoff = SignalEstimator(PolynomialDNI(3, 800u"ns"))
e = pickoff.(TrapezoidalChargeFilter(4.0u"µs", 2.0u"µs").(wvfs_pz), t50 .+ 5.0u"µs") ./ gain;

current_window = 15u"µs" .. 45u"µs"
wvfs_sg = SavitzkyGolayFilter(180u"ns", 3, 1).(wvfs_pz);

a_sg = get_a(wvfs_sg, current_window)
keys(a_sg[1])

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

wvfs_deriv = DerivativeFilter(1).(wvfs_pz)
wvfs_wsg = DerivativeFilter(1).(WeightedSavitzkyGolayFilter(length = 300u"ns", degree = 3, weightType = 2).(wvfs_pz))
wvfs_mwa = MovingWindowMultiFilter(100u"ns"; passes = 3).(wvfs_deriv)

a_wsg = get_a(wvfs_wsg, current_window)
a_mwa = get_a(wvfs_mwa, current_window)
a_raw = get_a(wvfs_deriv, current_window);

dep = findall(is_dep)
upper_width(x) = 100 * mad(x[x .> median(x)] .- median(x), center = 0, normalize = true) / median(x)

widths = map((raw = a_raw, sg = a_sg, wsg = a_wsg, mwa = a_mwa)) do a
    aoe = a.max[dep] ./ e[dep]
    upper_width(aoe)
end

aoe = a_sg.max ./ e
aoe ./= median(aoe[dep])

fig = Figure(size = (800, 380))
ax = Axis(fig[1, 1], xlabel = "A/E (normalized)", ylabel = "Events / 0.01", yscale = log10, limits = (0.4, 1.2, 0.5, nothing))
stephist!(ax, aoe, bins = 0.4:0.01:1.2, label = "all events")
stephist!(ax, aoe[dep], bins = 0.4:0.01:1.2, label = "DEP")
axislegend(ax, position = :lt)
fig

fig = Figure(size = (800, 450))
ax = Axis(fig[1, 1], xlabel = "Energy (keV)", ylabel = "A/E (normalized)", limits = (400, 2000, 0.4, 1.2))
h2 = fit(Histogram, (e, aoe), (400:40:2000, 0.4:0.02:1.2))
counts = replace(h2.weights, 0 => NaN)
hm = heatmap!(ax, h2.edges[1], h2.edges[2], counts, colorscale = log10)
Colorbar(fig[1, 2], hm, label = "Events")
fig

cut = 0.9
(survival_dep = mean(aoe[dep] .> cut), survival_continuum = mean(aoe[.!is_dep] .> cut))

t0 = get_t0(wvfs_pz, 80; mintot = 600u"ns")
qdrift = get_qdrift(wvfs_pz, t0, 2.5u"µs":0.1u"µs":5.0u"µs");

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
