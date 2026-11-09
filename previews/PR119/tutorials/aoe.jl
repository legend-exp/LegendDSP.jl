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

bl = signalstats.(wvfs, 0u"µs", 40u"µs")
wvfs_pz = InvCRFilter(τ).(shift_waveform.(wvfs, -bl.mean));

e_101010 = maximum.(TrapezoidalChargeFilter(10u"µs", 10u"µs").(wvfs_pz).signal)
t50 = get_threshold(wvfs_pz, 0.5 .* e_101010; mintot = 600u"ns")
pickoff = SignalEstimator(PolynomialDNI(3, 800u"ns"))
e = pickoff.(TrapezoidalChargeFilter(4.0u"µs", 2.0u"µs").(wvfs_pz), t50 .+ 5.0u"µs") ./ gain;

current_window = 40u"µs" .. 70u"µs"
wvfs_sg = SavitzkyGolayFilter(180u"ns", 3, 1).(wvfs_pz);

a_sg = get_a(wvfs_sg, current_window)
keys(a_sg[1])

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

wvfs_deriv = DerivativeFilter(1).(wvfs_pz)
wvfs_wsg = DerivativeFilter(1).(WeightedSavitzkyGolayFilter(length = 300u"ns", degree = 3, weightType = 2).(wvfs_pz))
wvfs_mwa = MovingWindowMultiFilter(100u"ns"; passes = 3).(wvfs_deriv)

a_wsg = get_a(wvfs_wsg, current_window)
a_mwa = get_a(wvfs_mwa, current_window)
a_raw = get_a(wvfs_deriv, current_window);

dep = findall(abs.(e .- 1592.5) .< 5)
upper_width(x) = 100 * mad(x[x .> median(x)] .- median(x), center = 0, normalize = true) / median(x)

widths = map((raw = a_raw, sg = a_sg, wsg = a_wsg, mwa = a_mwa)) do a
    aoe = a.max[dep] ./ e[dep]
    upper_width(aoe)
end

aoe = a_sg.max ./ e
aoe ./= median(aoe[dep])

bins = 0.4:0.01:1.2
fig = Figure(size = (800, 380))
ax = Axis(fig[1, 1], xlabel = "A/E (normalized)", ylabel = "Events / 0.01", yscale = log10, limits = (0.4, 1.2, 0.8, nothing))
for (selection, label) in ((eachindex(aoe), "all events"), (dep, "DEP window"))
    counts = fit(Histogram, aoe[selection], bins).weights
    stairs!(ax, bins, [max.(counts, 0.1); 0.1], step = :post, label = label)
end
axislegend(ax, position = :lt)
fig

fig = Figure(size = (800, 450))
ax = Axis(fig[1, 1], xlabel = "Energy (keV)", ylabel = "A/E (normalized)", limits = (400, 2700, 0.4, 1.2))
h2 = fit(Histogram, (e, aoe), (400:25:2700, 0.4:0.02:1.2))
hm = heatmap!(ax, h2.edges[1], h2.edges[2], replace(h2.weights, 0 => NaN), colorscale = log10)
Colorbar(fig[1, 2], hm, label = "Events")
fig

cut = 0.9
survival(lo, hi) = mean(aoe[lo .< e .< hi] .> cut)
(dep = survival(1588.5, 1596.5), sep = survival(2099.5, 2107.5), fep = survival(2610.5, 2618.5), continuum = survival(1700, 2050))

t0 = get_t0(wvfs_pz, 80; mintot = 600u"ns")
qdrift = get_qdrift(wvfs_pz, t0, 2.5u"µs":0.1u"µs":5.0u"µs");

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

n_lq = 600
e_lq_true = 1000 .+ 1600 .* rand(rng, n_lq)
late = [rand(rng) < 0.3 ? 0.02 + 0.06 * rand(rng) : 0.0 for _ in 1:n_lq]
wvfs_lq = simulate_aoe_waveforms(rng, e_lq_true, zeros(n_lq), fill(0.0u"ns", n_lq), late; gain = Float32(gain), τ = τ)

bl_lq = signalstats.(wvfs_lq, 0u"µs", 40u"µs")
wvfs_lq_pz = InvCRFilter(τ).(shift_waveform.(wvfs_lq, -bl_lq.mean))
e_lq = maximum.(TrapezoidalChargeFilter(10u"µs", 10u"µs").(wvfs_lq_pz).signal)
t80_lq = get_threshold(wvfs_lq_pz, 0.8 .* e_lq; mintot = 600u"ns");

early, gap, late_window = 2.5u"µs", 0.0u"µs", 2.5u"µs"
lq = get_lq(IntegratorFilter(1).(wvfs_lq_pz), t80_lq, early, gap, late_window);

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

lq_time = uconvert.(u"ns", lq .* step_t ./ e_lq)
has_late = late .> 0

fig = Figure(size = (800, 380))
ax = Axis(fig[1, 1], xlabel = "Energy (keV)", ylabel = "LQ / E (ns)", limits = (nothing, nothing, 0, 110))
scatter!(ax, e_lq[.!has_late] ./ gain, ustrip.(u"ns", lq_time[.!has_late]), markersize = 5, label = "no late charge")
scatter!(ax, e_lq[has_late] ./ gain, ustrip.(u"ns", lq_time[has_late]), markersize = 5, label = "late charge")
axislegend(ax, position = :lt)
fig
