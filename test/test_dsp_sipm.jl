using Test
using LegendDSP
using RadiationDetectorSignals
using RadiationDetectorDSP
using Unitful
using TypedTables
using PropDicts
using StructArrays
using ArraysOfArrays
using Random

function make_sipm_waveform(n=6250)
    t = (0:n-1) .* 16.0u"ns"
    signal = zeros(Float64, n)
    # flat baseline with a small SiPM-like pulse around 50µs
    pulse_start = round(Int, 50e3 / 16)  # ~3125
    pulse_width = 10
    amplitude = 5.0
    τ_samples = 30.0
    for i in 1:n
        if pulse_start <= i < pulse_start + pulse_width
            signal[i] = amplitude * (1 - exp(-(i - pulse_start) / 3.0))
        elseif i >= pulse_start + pulse_width
            signal[i] = amplitude * exp(-(i - pulse_start - pulse_width) / τ_samples)
        end
    end
    RDWaveform(t, signal)
end

function make_sipm_data(N=10)
    Table(
        waveform    = StructArray([make_sipm_waveform() for _ in 1:N]),
        baseline    = fill(0.0f0, N),
        timestamp   = fill(UInt64(0), N),
        eventnumber = UInt32.(1:N),
        daqenergy   = fill(UInt16(0), N),
    )
end

function make_sipm_config()
    PropDict(
        :t0_hpge_window => [47.0u"µs", 53.0u"µs"],
        :sg_flt_degree  => 3,
        :filters => PropDict(
            :sg => PropDict(
                :n_σ_threshold     => 3.0,
                :min_threshold     => -1.0,
                :max_threshold     => 1.0,
                :n_σ_dc_threshold  => 5.0,
                :min_dc_threshold  => -4.0,
                :max_dc_threshold  => 4.0,
                :min_tot_intersect => 70.0u"ns",
                :max_tot_intersect => 150.0u"ns",
            ),
            :trap => PropDict(
                :rt                => 100.0u"ns",
                :ft                => 50.0u"ns",
                :pz_tau            => 3.0u"µs",
                :n_σ_threshold     => 3.5,
                :min_threshold     => -1.5,
                :max_threshold     => 1.5,
                :n_σ_dc_threshold  => 5.0,
                :min_dc_threshold  => -3.0,
                :max_dc_threshold  => 3.0,
                :min_tot_intersect => 48.0u"ns",
                :max_tot_intersect => 250.0u"ns",
            ),
        ),
    )
end

@testset "dsp_sipm" begin
    data = make_sipm_data()
    config = make_sipm_config()
    pars_opt = PropDict(:sg => PropDict(:wl => 200.0u"ns"))

    result = dsp_sipm(data, config, pars_opt)

    # check that result is a Table with the right number of rows
    @test result isa TypedTables.Table
    @test length(result) == 10

    # check all expected output columns exist
    expected_keys = [:blfc, :timestamp, :eventID_fadc, :e_fc,
        :t_max, :t_min, :t_max_lar, :t_min_lar,
        :e_max, :e_min, :e_max_lar, :e_min_lar,
        :blmean, :blsigma, :blslope, :bloffset,
        :wfmean, :wfsigma, :wfslope, :wfoffset,
        :threshold_sg, :threshold_DC_sg,
        :trig_pos_sg, :trig_max_sg, :trig_pos_DC_sg, :trig_max_DC_sg,
        :threshold_trap, :threshold_DC_trap,
        :trig_pos_trap, :trig_pos_high_trap, :trig_pos_tot_trap, :trig_max_trap,
        :trig_pos_DC_trap, :trig_pos_high_DC_trap, :trig_pos_tot_DC_trap, :trig_max_DC_trap]
    for k in expected_keys
        @test hasproperty(result, k)
    end

    # basic sanity: timestamps and event IDs passed through
    @test all(result.timestamp .== 0)
    @test result.eventID_fadc == UInt32.(1:10)

    # thresholds should be finite and non-negative
    @test all(isfinite.(result.threshold_sg))
    @test all(result.threshold_sg .>= 0)
    @test all(isfinite.(result.threshold_trap))
    @test all(result.threshold_trap .>= 0)

    # time values should be in µs and within waveform range
    @test all(0.0u"µs" .<= result.t_max .<= 100.0u"µs")
    @test all(0.0u"µs" .<= result.t_min .<= 100.0u"µs")
end

# Matrix-backed signals with a Float32 time axis, as returned by decode_data,
# and a mix of events with and without triggers.
function make_mixed_sipm_data(N=10)
    rng = MersenneTwister(1)
    pulse = make_sipm_waveform().signal
    t = range(0f0u"ns", step=16f0u"ns", length=length(pulse))
    signals = reduce(hcat, [(isodd(i) ? pulse : zero(pulse)) .+ 0.01 .* randn(rng, length(pulse)) for i in 1:N])
    Table(
        waveform    = ArrayOfRDWaveforms((fill(t, N), VectorOfSimilarVectors(signals))),
        baseline    = fill(0.0f0, N),
        timestamp   = fill(UInt64(0), N),
        eventnumber = UInt32.(1:N),
        daqenergy   = fill(UInt16(0), N),
    )
end

@testset "dsp_sipm baseline window" begin
    data = make_mixed_sipm_data()
    config = make_sipm_config()
    result = dsp_sipm(data, config, PropDict(:sg => PropDict(:wl => 200.0u"ns")))

    has_trig = .!isempty.(result.trig_pos_sg)
    @test any(has_trig) && !all(has_trig)

    # without a trigger, the baseline window is the full waveform
    @test result.blmean[.!has_trig] ≈ result.wfmean[.!has_trig]
    @test result.blsigma[.!has_trig] ≈ result.wfsigma[.!has_trig]

    # with a trigger, the baseline window ends at the first trigger
    wvfs = LegendDSP.shift_waveform.(data.waveform, 0.0)
    wvfs_int = IntegratorFilter(gain=1).(SavitzkyGolayFilter(200.0u"ns", config.sg_flt_degree, 1).(wvfs))
    t_min = first(wvfs_int[1].time)
    for i in findall(has_trig)
        @test result.blmean[i] ≈ signalstats(wvfs_int[i], t_min, minimum(result.trig_pos_sg[i])).mean
    end
end

# Note: dsp_sipm_compressed is not tested here because it requires
# bit-drop compressed waveform data (waveform_bit_drop) decoded via
# LegendDataTypes.decode_data, which is hard to mock with synthetic data.
# The DSP logic is identical to dsp_sipm.
