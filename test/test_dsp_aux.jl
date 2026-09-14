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
        baseline = Float32[100, 40, -8],
        timestamp = UInt64[11, 12, 13],
        eventnumber = UInt32[21, 22, 23],
        daqenergy = UInt16[31, 32, 33],
    )
end

@testset "auxiliary DSP" begin
    data = make_aux_data()
    config = make_fake_config()
    expected_e_max = [125, 72, -5]
    expected_t_max = [32, 16, 16]u"ns"

    @testset "uncompressed waveforms" begin
        result = aux_dsp(data, config)

        @test result isa TypedTables.Table
        @test length(result) == 3
        @test columnnames(result) == (
            :e_max, :t_max, :blfc, :timestamp, :eventID_fadc, :e_fc,
        )
        @test result.e_max == expected_e_max
        @test result.t_max == expected_t_max
        @test result.blfc == data.baseline
        @test result.timestamp == data.timestamp
        @test result.eventID_fadc == data.eventnumber
        @test result.e_fc == data.daqenergy
    end

    @testset "compressed waveforms" begin
        result = aux_dsp_compressed(data, config)

        @test result isa TypedTables.Table
        @test length(result) == 3
        @test columnnames(result) == (
            :e_max, :t_max, :blfc, :timestamp, :eventID_fadc, :e_fc,
        )
        @test result.e_max == expected_e_max
        @test result.t_max == expected_t_max
        @test result.blfc == data.baseline
        @test result.timestamp == data.timestamp
        @test result.eventID_fadc == data.eventnumber
        @test result.e_fc == data.daqenergy
    end

    @testset "maxima use raw samples" begin
        # A baseline subtraction would produce 25, 32, and 3 instead.
        @test aux_dsp(data, config).e_max == expected_e_max
        @test aux_dsp_compressed(data, config).e_max == expected_e_max
    end
end

