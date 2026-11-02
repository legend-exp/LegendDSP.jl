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
            :mwa => PropDict(
                :wl_trigger        => 160.0u"ns",
                :passes            => 4,
                :n_σ_threshold     => 3.0,
                :min_threshold     => -1.0,
                :max_threshold     => 1.0,
                :solve_overlap     => true,
                :merge_distance    => 0.25,
            ),
        ),
    )
end

make_sipm_pars() = PropDict(:sg => PropDict(:wl => 200.0u"ns"), :mwa => PropDict(:wl => 800.0u"ns"))

# one photon = an exponential pulse of height A and decay ratio r per sample; its exact inverse kernel is
# y[i] = (x[i] - r x[i-1]) / A, i.e. the taps (lag -1, 0, +1) = (0, 1/A, -r/A), which turns it into a unit stick
const PHOTON_A, PHOTON_R = 8.0, exp(-1 / 40)
photon_pulse(n, i0) = [i < i0 ? 0.0 : PHOTON_A * PHOTON_R^(i - i0) for i in 1:n]
make_sipm_kernel_pars() = PropDict(:kernel => PropDict(:dt => 16.0u"ns", :taps => [0.0, 1 / PHOTON_A, -PHOTON_R / PHOTON_A]))

# photons at 50 µs (one), 52 µs (two, 300 ns apart) on the odd rows, noise only on the even rows
function make_photon_data(N = 10)
    rng = MersenneTwister(2)
    n = 6250
    t = range(0f0u"ns", step = 16f0u"ns", length = n)
    i1 = round(Int, 50e3 / 16)
    i2 = round(Int, 52e3 / 16)
    signals = reduce(hcat, [(isodd(k) ? photon_pulse(n, i1) .+ photon_pulse(n, i2) .+ photon_pulse(n, i2 + 19) : zeros(n)) .+ 0.05 .* randn(rng, n) for k in 1:N])
    Table(
        waveform    = ArrayOfRDWaveforms((fill(t, N), VectorOfSimilarVectors(signals))),
        baseline    = fill(0.0f0, N),
        timestamp   = fill(UInt64(0), N),
        eventnumber = UInt32.(1:N),
        daqenergy   = fill(UInt16(0), N),
    )
end

@testset "dsp_sipm" begin
    data = make_sipm_data()
    config = make_sipm_config()
    result = dsp_sipm(data, config, make_sipm_pars())

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
        :trig_pos_DC_trap, :trig_pos_high_DC_trap, :trig_pos_tot_DC_trap, :trig_max_DC_trap,
        :threshold_mwa, :trig_pos_mwa, :trig_pos_tot_mwa, :trig_max_mwa]
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
    result = dsp_sipm(data, config, make_sipm_pars())

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

@testset "dsp_sipm MWA chain" begin
    data = make_photon_data()
    config = make_sipm_config()
    has_pulse = isodd.(1:length(data))
    pars = merge(make_sipm_pars(), make_sipm_kernel_pars())

    # matrix method: sites on the 160-ns trace, amplitudes solved on the 800-ns trace; the two photons 300 ns apart
    # overlap on the 800-ns trace and are separated by the solve
    result = dsp_sipm(data, config, pars)
    for i in findall(has_pulse)
        big = findall(>(0.5), result.trig_max_mwa[i])
        @test length(big) == 3
        @test result.trig_max_mwa[i][big] ≈ [1.0, 1.0, 1.0] atol = 0.1
        @test result.trig_pos_mwa[i][big[1]] ≈ 50.0u"µs" atol = 0.1u"µs"
        @test all(result.trig_pos_tot_mwa[i][big] .> 0.0u"ns")
    end
    for i in findall(.!has_pulse)
        @test all(abs.(result.trig_max_mwa[i]) .< 0.3)
    end
    @test eltype(eltype(result.trig_pos_mwa)) == eltype(first(data.waveform).time)

    # without the matrix method the sites are the extrema of the 800-ns trace: the close pair is one site of about two photons
    config_nosolve = deepcopy(config)
    config_nosolve.filters.mwa.solve_overlap = false
    result_ns = dsp_sipm(data, config_nosolve, pars)
    for i in findall(has_pulse)
        big = findall(>(0.5), result_ns.trig_max_mwa[i])
        @test length(big) == 2
        @test result_ns.trig_max_mwa[i][big[1]] ≈ 1.0 atol = 0.1
        @test result_ns.trig_max_mwa[i][big[2]] > 1.5
    end

    # without template pars the derivative takes the place of the kernel: same columns, amplitudes in ADC
    result_nk = dsp_sipm(data, config, make_sipm_pars())
    @test propertynames(result_nk) == propertynames(result)
    @test all(!isempty(result_nk.trig_pos_mwa[i]) for i in findall(has_pulse))

    # the MWA chain needs an even number of passes and the kernel sampling must match the waveforms
    config_odd = deepcopy(config)
    config_odd.filters.mwa.passes = 3
    @test_throws ArgumentError dsp_sipm(data, config_odd, pars)
    pars_dt = deepcopy(pars)
    pars_dt.kernel.dt = 8.0u"ns"
    @test_throws ArgumentError dsp_sipm(data, config, pars_dt)
end
