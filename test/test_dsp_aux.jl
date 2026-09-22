using Test
using LegendDSP
using RadiationDetectorSignals
using StructArrays
using TypedTables
using Unitful

function make_aux_waveform(signal, step_size = 16u"ns")
    RDWaveform(range(0u"ns", step = step_size, length = length(signal)), signal)
end

function make_aux_data()
    waveforms = [
        make_aux_waveform([100, 101, 125, 103, 102]),
        make_aux_waveform([40, 72, 45, 44, 43]),
        make_aux_waveform([-8, -5, -6, -7, -9]),
    ]

    Table(
        waveform = StructArray(waveforms),
        waveform_presummed = StructArray(waveforms),
        presum_rate = UInt8[8, 4, 2],
        baseline = Float32[100, 40, -8],
        timestamp = UInt64[11, 12, 13],
        eventnumber = UInt32[21, 22, 23],
        daqenergy = UInt16[31, 32, 33],
    )
end

@testset "auxiliary DSP" begin
    data = make_aux_data()
    config = make_fake_config()

    # Decoded maxima are [125, 72, -5]. Uncompressed data subtracts the
    # baseline directly, while compressed data first scales it by presum_rate.
    expected_uncompressed_e_max = Float32[25, 32, 3]
    expected_compressed_e_max = Float32[-675, -88, 11]
    expected_t_max = [32, 16, 16]u"ns"

    @testset "uncompressed waveforms" begin
        result = dsp_aux(data, config)

        @test result isa TypedTables.Table
        @test length(result) == 3
        @test columnnames(result) == (
            :e_max, :t_max, :blfc, :timestamp, :eventID_fadc, :e_fc,
        )
        @test result.e_max == expected_uncompressed_e_max
        @test result.t_max == expected_t_max
        @test result.blfc == data.baseline
        @test result.timestamp == data.timestamp
        @test result.eventID_fadc == data.eventnumber
        @test result.e_fc == data.daqenergy
    end

    @testset "compressed waveforms" begin
        result = dsp_aux_compressed(data, config)

        @test result isa TypedTables.Table
        @test length(result) == 3
        @test columnnames(result) == (
            :e_max, :t_max, :blfc, :timestamp, :eventID_fadc, :e_fc,
        )
        @test result.e_max == expected_compressed_e_max
        @test result.t_max == expected_t_max
        @test result.blfc == data.baseline
        @test result.timestamp == data.timestamp
        @test result.eventID_fadc == data.eventnumber
        @test result.e_fc == data.daqenergy
    end

    @testset "maxima include baseline subtraction" begin
        @test dsp_aux(data, config).e_max == expected_uncompressed_e_max
        @test dsp_aux_compressed(data, config).e_max == expected_compressed_e_max
    end
end

@testset "pulser DSP" begin
    compressed_data = make_fake_data(3)
    data = Table(
        waveform = compressed_data.waveform_presummed,
        baseline = compressed_data.baseline,
        timestamp = compressed_data.timestamp,
        eventnumber = compressed_data.eventnumber,
        daqenergy = compressed_data.daqenergy,
    )
    config = make_fake_config()

    for (name, result) in (
        "uncompressed waveforms" => dsp_puls(data, config),
        "compressed waveforms" => dsp_puls_compressed(compressed_data, config),
    )
        @testset "$name" begin
            @test result isa TypedTables.Table
            @test length(result) == 3
            @test columnnames(result) == (
                :blmean, :blsigma, :blslope, :bloffset, :t50,
                :e_max, :e_10410, :blfc, :timestamp, :eventID_fadc, :e_fc,
            )
            @test all(isapprox.(result.blmean, 1000; atol = 1e-10))
            @test all(isfinite, result.blsigma)
            @test all(isfinite, result.blslope)
            @test all(isfinite, result.bloffset)
            @test all(48u"μs" .< result.t50 .< 50u"μs")
            @test all(result.e_max .> 9_900)
            @test all(result.e_10410 .> 0)
            @test result.blfc == data.baseline
            @test result.timestamp == data.timestamp
            @test result.eventID_fadc == data.eventnumber
            @test result.e_fc == data.daqenergy
        end
    end
end
