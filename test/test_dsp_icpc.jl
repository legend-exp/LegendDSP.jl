using Test
using LegendDSP
using RadiationDetectorSignals
using RadiationDetectorDSP
using Unitful
using TypedTables
using PropDicts
using IntervalSets
using StructArrays

function make_fake_waveform(n=8192, presum_rate=1)
    t = (0:n-1) .* 16.0u"ns"
    signal = zeros(Float64, n)
    
    baseline_end = round(Int, 48e3 / 16)
    rise_end     = round(Int, 50e3 / 16)
    amplitude    = 10000.0 * presum_rate
    τ_samples    = 500e3 / 16
    baseline_offset = 1000.0 

    for i in 1:n
        if i < baseline_end
            signal[i] = baseline_offset
        elseif i < rise_end
            signal[i] = baseline_offset + amplitude * (i - baseline_end) / (rise_end - baseline_end)
        else
            signal[i] = baseline_offset + amplitude * exp(-(i - rise_end) / τ_samples)
        end
    end

    RDWaveform(t, signal)
end


function make_fake_data(N=3)
    Table(
        waveform_presummed = StructArray([make_fake_waveform() for _ in 1:N]),
        waveform_windowed  = StructArray([make_fake_waveform() for _ in 1:N]),
        presum_rate        = fill(UInt16(1), N),
        baseline           = fill(0.0f0, N),
        timestamp          = fill(UInt64(0), N),
        eventnumber        = UInt32.(1:N),
        daqenergy          = fill(UInt16(0), N),
        t_sat_lo           = fill(UInt16(0), N),
        t_sat_hi           = fill(UInt16(0), N),
        deadtime           = fill(UInt16(0), N),
    )
end

function make_fake_config()
    pd = PropDict(

        :enc_pickoff_trap => 40.0u"µs",
        :enc_pickoff_zac  => 41.0u"µs",
        :enc_pickoff_cusp => 41.0u"µs",

        :bl_window => PropDict(
            :min => 0.0u"µs",
            :max => 39.0u"µs",
        ),
        :tail_window => PropDict(
            :min => 70.0u"µs",
            :max => 110.0u"µs",
        ),
        :current_window => PropDict(
            :min => 43.0u"µs",
            :max => 62.0u"µs",
        ),

        :auxbl1_window => PropDict(
            :min => 0.0u"µs",
            :max => 20.0u"µs",
        ),
        :auxbl2_window => PropDict(
            :min => 20.0u"µs",
            :max => 39.0u"µs",
        ),

        :auxpz1_window => PropDict(
            :min => 70.0u"µs",
            :max => 90.0u"µs",
        ),
        :auxpz2_window => PropDict(
            :min => 90.0u"µs",
            :max => 110.0u"µs",
        ),
        :flt_length_cusp => 38.0u"µs",
        :flt_length_zac  => 38.0u"µs",
        :t0_threshold             => 4.0,
        :inTraceCut_std_threshold => 5,
        :sg_flt_degree            => 3,
        :qdrift_int_length => (2.5u"µs":0.1u"µs":5.0u"µs"),
        :lq_int_length     => (2.5u"µs":0.1u"µs":5.0u"µs"),

        :e_grid_trap => PropDict(
            :rt => PropDict(
                :start => 1.0u"µs",
                :stop  => 16.0u"µs",
                :step  => 0.5u"µs",
            ),
            :ft => PropDict(
                :start => 1.0u"µs",
                :stop  => 4.0u"µs",
                :step  => 0.2u"µs",
            ),
        ),

        :e_grid_zac => PropDict(
            :rt => PropDict(
                :start => 1.0u"µs",
                :stop  => 16.0u"µs",
                :step  => 0.5u"µs",
            ),
            :ft => PropDict(
                :start => 1.0u"µs",
                :stop  => 4.0u"µs",
                :step  => 0.2u"µs",
            ),
        ),

        :e_grid_cusp => PropDict(
            :rt => PropDict(
                :start => 1.0u"µs",
                :stop  => 16.0u"µs",
                :step  => 0.5u"µs",
            ),
            :ft => PropDict(
                :start => 1.0u"µs",
                :stop  => 4.0u"µs",
                :step  => 0.2u"µs",
            ),
        ),

        :a_grid_wl_sg => PropDict(
            :start => 30.0u"ns",
            :stop  => 350.0u"ns",
            :step  => 32.0u"ns",
        ),

        :flt_defaults => PropDict(
            :sg   => 100.0u"ns",
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


@testset "dsp_icpc_compressed" begin
    data   = make_fake_data(3)
    config = make_fake_config()
    τ      = 500u"µs"
    pars_filter = PropDict()  # falls back to flt_defaults

    result = dsp_icpc_compressed(data, config, τ, pars_filter)

    @testset "Output table shape" begin
        @test result isa TypedTables.Table
        @test length(result) == 3
        # Spot-check a few key columns are present
        cols = columnnames(result)
        for col in [:blmean, :blsigma, :blslope, :bloffset,
                    :tailmean, :tailsigma, :tailslope, :tailoffset,
                    :t0, :t50, :t90, :drift_time,
                    :e_10410, :e_313, :e_trap, :e_cusp, :e_zac,
                    :qdrift, :lq, :a_sg,
                    :n_sat_low, :n_sat_high,
                    :inTrace_intersect, :inTrace_n,
                    :e_10410_inv, :e_313_inv, :t0_inv]
            @test col in cols
        end
    end

    @testset "Timing ordering" begin
        @test all(result.t0  .< result.t50)
        @test all(result.t50 .< result.t90)
        @test all(ustrip.(result.drift_time) .>= 0)
    end

    @testset "Energies finite" begin
        for col in [:e_10410, :e_313, :e_trap]
            @test all(isfinite, ustrip.(getproperty(result, col)))
        end
    end
end
