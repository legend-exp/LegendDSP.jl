# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

using Test
using LegendDSP
using RadiationDetectorSignals
using Statistics
using Unitful

@testset "Test thresholdstats" begin

    # thresholdstats is supposed to do the following, but faster:
    compare(wf::RDWaveform{<:Any,T}, _min = -Inf*unit(T), _max = Inf*unit(T)) where {T} = std(wf.signal[_min .<= wf.signal .<= _max])
    compare(wf::AbstractVector{T}, _min = -Inf*unit(T), _max = Inf*unit(T)) where {T} = std(wf[_min .<= wf .<= _max])

    @testset "Test single waveform" begin
        @testset "No units" begin
            σ = 10.0 * rand()
            wf = RDWaveform((1:1000)u"ns", σ .* randn(10000))
            # test that the implementation returns the initial σ within 5%
            @test isapprox(thresholdstats(wf), σ, rtol = 0.05)
            # test that the result agrees with the "slow" implementation within 0.1%
            @test isapprox(thresholdstats(wf), compare(wf), rtol = 0.001)
            # test that the implementation of min and max work as intended
            @test(all([begin
                _min = -σ * rand(); _max =  σ * rand()
                isapprox(thresholdstats(wf, _min, _max), compare(wf, _min, _max), rtol = 0.005)
                end for _ in 1:1000
            ]))
        end
        @testset "With units" begin
            σ = 10.0u"pC" * rand()
            wf = RDWaveform((1:1000)u"ns", σ .* randn(10000))
            # test that the implementation returns the initial σ within 5%
            @test isapprox(thresholdstats(wf), σ, rtol = 0.05)
            # test that the result agrees with the "slow" implementation within 0.1%
            @test isapprox(thresholdstats(wf), compare(wf), rtol = 0.001)
            # test that the implementation of min and max work as intended
            @test(all([begin
                _min = -σ * rand(); _max =  σ * rand()
                isapprox(thresholdstats(wf, _min, _max), compare(wf, _min, _max), rtol = 0.005)
                end for _ in 1:1000
            ]))
        end
    end

    @testset "Test single array" begin
        @testset "No units" begin
            σ = 10.0 * rand()
            wf = σ .* randn(10000)
            # test that the implementation returns the initial σ within 5%
            @test isapprox(thresholdstats(wf), σ, rtol = 0.05)
            # test that the result agrees with the "slow" implementation within 0.1%
            @test isapprox(thresholdstats(wf), compare(wf), rtol = 0.001)
            # test that the implementation of min and max work as intended
            @test(all([begin
                _min = -σ * rand(); _max =  σ * rand()
                isapprox(thresholdstats(wf, _min, _max), compare(wf, _min, _max), rtol = 0.005)
                end for _ in 1:1000
            ]))
        end
    end
end
