# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

using Test
using LegendDSP
using RadiationDetectorSignals
using RadiationDetectorDSP
using Unitful
using PropDicts
using IntervalSets
using Random

# Grids are small to keep the test fast; the pulse rises linearly at
# 48 to 52 µs (jittered per event), decays with τ = 500 µs and carries
# Gaussian noise. The baseline window ends before the earliest onset.
function fltopt_config()
    pd = PropDict(
        :enc_pickoff_trap => 40.0u"µs",
        :enc_pickoff_zac  => 41.0u"µs",
        :enc_pickoff_cusp => 41.0u"µs",
        :bl_window      => PropDict(:min => 0.0u"µs",  :max => 39.0u"µs"),
        :tail_window    => PropDict(:min => 70.0u"µs", :max => 110.0u"µs"),
        :current_window => PropDict(:min => 43.0u"µs", :max => 62.0u"µs"),
        :auxpz1_window  => PropDict(:min => 70.0u"µs", :max => 90.0u"µs"),
        :auxpz2_window  => PropDict(:min => 90.0u"µs", :max => 110.0u"µs"),
        :flt_length_cusp => 38.0u"µs",
        :flt_length_zac  => 38.0u"µs",
        :t0_threshold             => 4.0,
        :inTraceCut_std_threshold => 5,
        :sg_flt_degree            => 3,
        :qdrift_int_length => (2.5u"µs":0.1u"µs":5.0u"µs"),
        :lq_int_length     => (2.5u"µs":0.1u"µs":5.0u"µs"),
        :e_grid_trap => PropDict(
            :rt => PropDict(:start => 2.0u"µs", :stop => 6.0u"µs", :step => 2.0u"µs"),
            :ft => PropDict(:start => 1.0u"µs", :stop => 3.0u"µs", :step => 0.5u"µs"),
        ),
        :e_grid_zac => PropDict(
            :rt => PropDict(:start => 2.0u"µs", :stop => 6.0u"µs", :step => 2.0u"µs"),
            :ft => PropDict(:start => 1.0u"µs", :stop => 3.0u"µs", :step => 0.5u"µs"),
        ),
        :e_grid_cusp => PropDict(
            :rt => PropDict(:start => 2.0u"µs", :stop => 6.0u"µs", :step => 2.0u"µs"),
            :ft => PropDict(:start => 1.0u"µs", :stop => 3.0u"µs", :step => 0.5u"µs"),
        ),
        :a_grid_wl_sg => PropDict(:start => 30.0u"ns", :stop => 350.0u"ns", :step => 32.0u"ns"),
        :a_grid_wl_mwa => PropDict(:start => 16.0u"ns", :stop => 576.0u"ns", :step => 16.0u"ns"),
        :flt_defaults => PropDict(
            :sg   => 100.0u"ns",
            :mwa  => 96.0u"ns",
            :trap => PropDict(:rt => 5.0u"µs", :ft => 2.5u"µs"),
            :zac  => PropDict(:rt => 5.0u"µs", :ft => 2.5u"µs"),
            :cusp => PropDict(:rt => 5.0u"µs", :ft => 2.5u"µs"),
        ),
        :kwargs_pars => PropDict(
            :fc_bit_depth             => 16,
            :t0_flt_pars              => [40.0u"ns", 100.0u"ns", 2000.0u"ns"],
            :t0_mintot                => 1500.0u"ns",
            :tx_mintot                => 32.0u"ns",
            :intrace_mintot           => 100.0u"ns",
            :int_interpolation_order  => 3,
            :int_interpolation_length => 100.0u"ns",
            :sig_interpolation_order  => 3,
            :sig_interpolation_length => 700.0u"ns",
        ),
    )
    return DSPConfig(pd)
end

# Returns the waveform set and the per-event onset times in µs.
function fltopt_waveforms(rng, N; n=8192, dt=16.0u"ns", amplitude=10000.0,
                          τ=500.0u"µs", σ=5.0, baseline=1000.0, rise=2.0)
    t = (0:n-1) .* dt
    ts = ustrip.(u"µs", t)
    τ_µs = ustrip(u"µs", τ)
    onsets = 48.0 .+ 4.0 .* rand(rng, N)
    wvfs = map(onsets) do t_on
        signal = map(ts) do x
            x < t_on          ? baseline :
            x < t_on + rise   ? baseline + amplitude * (x - t_on) / rise :
                                baseline + amplitude * exp(-(x - t_on - rise) / τ_µs)
        end
        RDWaveform(t, signal .+ σ .* randn(rng, n))
    end
    return ArrayOfRDWaveforms(wvfs), onsets
end

@testset "single-point filter optimization" begin
    config = fltopt_config()
    τ = 500.0u"µs"
    wvfs, onsets = fltopt_waveforms(Xoshiro(20260929), 50)
    prep = dsp_flt_optimization_prepare(wvfs, config, τ)
    rt = 4.0u"µs"
    grids = (
        trap = (dsp_trap_ft_optimization, config.e_grid_ft_trap),
        cusp = (dsp_cusp_ft_optimization, config.e_grid_ft_cusp),
        zac  = (dsp_zac_ft_optimization,  config.e_grid_ft_zac),
    )

    @testset "prepare" begin
        @test keys(prep) == (:wvfs_pz, :t50, :bl_stats, :signal_estimator, :sampling_period)
        @test length(prep.wvfs_pz) == length(wvfs)
        @test prep.sampling_period == 16.0u"ns"
        # t50 is the 50% crossing of the linear rise, one µs after the onset
        @test all(abs.(ustrip.(u"µs", prep.t50) .- (onsets .+ 1.0)) .< 0.2)
        @test all(abs.(prep.bl_stats.mean .- 1000.0) .< 5.0)
    end

    @testset "dsp_energy_rtft equals the matching grid row" begin
        for flt in (:trap, :cusp, :zac)
            f_grid, ft_grid = grids[flt]
            e_grid = f_grid(wvfs, config, τ, rt)
            for ft_index in (firstindex(ft_grid), lastindex(ft_grid))
                e = dsp_energy_rtft(prep, config, rt, ft_grid[ft_index], flt)
                @test length(e) == length(wvfs)
                @test all(isfinite, e)
                @test isequal(Float32.(e), e_grid[ft_index, :])
            end
        end
    end

    @testset "explicit Val and Symbol dispatch agree" begin
        ft = 2.0u"µs"
        for flt in (:trap, :cusp, :zac)
            @test isequal(dsp_energy_rtft(prep, config, rt, ft, Val(flt)),
                          dsp_energy_rtft(prep, config, rt, ft, flt))
        end
    end

    @testset "unknown filter" begin
        @test_throws MethodError dsp_energy_rtft(prep, config, rt, 2.0u"µs", :unknown)
    end

    @testset "threaded evaluation equals serial evaluation" begin
        points = [(rt, ft, flt) for flt in (:trap, :cusp, :zac), rt in (2.0u"µs", 4.0u"µs"), ft in (1.0u"µs", 2.5u"µs")]
        serial = [dsp_energy_rtft(prep, config, p...) for p in points]
        threaded = similar(serial)
        Threads.@threads for i in eachindex(points)
            threaded[i] = dsp_energy_rtft(prep, config, points[i]...)
        end
        @test all(isequal.(serial, threaded))
    end

    @testset "grid functions equal the per-ft reference" begin
        # The reference is computed from the public filters only, independently
        # of the functions under test.
        function reference_grid(wvfs, config, τ, rt, ft_grid, make_flt, pickoff)
            bl_stats = signalstats.(wvfs, leftendpoint(config.bl_window), rightendpoint(config.bl_window))
            wvfs_pz = InvCRFilter(τ).(shift_waveform.(wvfs, -bl_stats.mean))
            signal_estimator = SignalEstimator(PolynomialDNI(config.kwargs_pars.sig_interpolation_order, config.kwargs_pars.sig_interpolation_length))
            t50 = get_threshold(wvfs_pz, maximum.(wvfs_pz.signal) .* 0.5; mintot=config.kwargs_pars.tx_mintot)
            e_grid = Array{Union{Missing, Float32}}(missing, length(ft_grid), length(wvfs_pz))
            for (f, ft) in pairs(ft_grid)
                e_grid[f, :] = signal_estimator.(make_flt(rt, ft).(wvfs_pz), pickoff(t50, rt, ft))
            end
            return e_grid
        end

        scale(flt_length) = ustrip(NoUnits, flt_length / step(first(wvfs).time))
        refs = (
            trap = ((rt, ft) -> TrapezoidalChargeFilter(rt, ft),
                    (t50, rt, ft) -> t50 .+ (rt + ft/2)),
            cusp = ((rt, ft) -> CUSPChargeFilter(rt, ft, 10000000.0u"µs", config.flt_length_cusp, scale(config.flt_length_cusp)),
                    (t50, rt, ft) -> t50 .+ config.flt_length_cusp/2),
            zac  = ((rt, ft) -> ZACChargeFilter(rt, ft, 10000000.0u"µs", config.flt_length_zac, scale(config.flt_length_zac)),
                    (t50, rt, ft) -> t50 .+ config.flt_length_zac/2),
        )
        for flt in (:trap, :cusp, :zac)
            f_grid, ft_grid = grids[flt]
            e_new = f_grid(wvfs, config, τ, rt)
            e_ref = reference_grid(wvfs, config, τ, rt, ft_grid, refs[flt]...)
            @test typeof(e_new) == typeof(e_ref)
            @test size(e_new) == (length(ft_grid), length(wvfs))
            @test isequal(e_new, e_ref)
        end
    end

    @testset "pulse shape" begin
        shape = dsp_pulse_shape_flt_optimization(prep, config, 100.0u"ns")
        @test keys(shape) == (:t10, :t90, :a)
        @test length(shape.t10) == length(shape.t90) == length(shape.a) == length(wvfs)
        @test all(shape.t90 .> shape.t10)
        @test all(isfinite, ustrip.(shape.a))
        # 10% and 90% crossings of the 2 µs linear rise lie 0.2 and 1.8 µs after the onset
        @test all(abs.(ustrip.(u"µs", shape.t10) .- (onsets .+ 0.2)) .< 0.2)
        @test all(abs.(ustrip.(u"µs", shape.t90) .- (onsets .+ 1.8)) .< 0.2)
        # the SG current maximum equals the direct computation on wvfs_pz
        a_ref = get_wvf_maximum.(SavitzkyGolayFilter(100.0u"ns", config.sg_flt_degree, 1).(prep.wvfs_pz),
                                 leftendpoint(config.current_window), rightendpoint(config.current_window)).max
        @test isequal(shape.a, a_ref)
    end
end
