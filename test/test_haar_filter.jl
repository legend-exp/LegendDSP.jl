using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful

@testset "Test HaarTypeWaveletFilter" begin

    T = Float64

    @testset "Step waveform" begin
        
        wvf = RDWaveform(range(0u"μs", step = 32u"μs", length = 256), vcat(ones(T,128), 2 .* ones(T,128)))
        
        @testset "HaarTypeWaveletFilter(2)" begin

            flt = HaarTypeWaveletFilter(2)

            # Apply HaarTypeWaveletFilter once
            wvf_out = flt(wvf)
            @test wvf_out.time == range(0u"μs", step = 64u"μs", length = 128)
            @test isapprox(wvf_out.signal, vcat(fill(T(sqrt(2)), 64), fill(T(sqrt(8)), 64)), rtol = eps(T))

            # Apply HaarTypeWaveletFilter multiple times
            wvf_out = flt(wvf_out)
            @test wvf_out.time == range(0u"μs", step = 128u"μs", length = 64)
            @test isapprox(wvf_out.signal, vcat(fill(T(2), 32), fill(T(4), 32)), rtol = eps(T))
        end

        @testset "HaarTypeWaveletFilter(4)" begin

            flt = HaarTypeWaveletFilter(4)

            # Apply HaarTypeWaveletFilter once
            wvf_out = flt(wvf)
            @test wvf_out.time == range(0u"μs", step = 128u"μs", length = 64)
            @test isapprox(wvf_out.signal, vcat(fill(T(sqrt(2)), 32), fill(T(sqrt(8)), 32)), rtol = eps(T))
        end
    end

    @testset "Linear waveform" begin

        wvf = RDWaveform(range(0u"μs", step = 32u"μs", length = 256), collect(T, 0:255))
    
        @testset "HaarTypeWaveletFilter(2)" begin

            flt = HaarTypeWaveletFilter(2)

            # Apply HaarTypeWaveletFilter once
            wvf_out = flt(wvf)
            @test wvf_out.time == range(0u"μs", step = 64u"μs", length = 128)
            @test isapprox(wvf_out.signal, collect(T, (0.5:2:254.5).*sqrt(2)), rtol = eps(T))

            # Apply HaarTypeWaveletFilter multiple times
            wvf_out = flt(wvf_out)
            @test wvf_out.time == range(0u"μs", step = 128u"μs", length = 64)
            @test isapprox(wvf_out.signal, collect(T, 3:8:510))
        end

        @testset "HaarTypeWaveletFilter(4)" begin

            flt = HaarTypeWaveletFilter(4)

            # Apply HaarTypeWaveletFilter once
            wvf_out = flt(wvf)
            @test wvf_out.time == range(0u"μs", step = 128u"μs", length = 64)
            @info diff(wvf_out.signal)
            @test isapprox(wvf_out.signal, collect(T, (0.5:4:252.5)*sqrt(2)), rtol = 10*eps(T))
        end
    end

    # TODO: add tests for some more realistic waveforms
     
end
