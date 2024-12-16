using Test
using LegendDSP
using RadiationDetectorSignals: RDWaveform
using Unitful

@testset "Test DerivativeFilter" begin
    T = Float64
    sig = rand(100)
    wvf = RDWaveform(range(0u"ns", step = 16u"ns", length = length(sig)), sig)

    flt = DerivativeFilter()
    wvf_out = flt(wvf)
    @test isapprox(wvf_out.signal, vcat(sig[2] - sig[1], diff(sig)...))

    gain = rand()
    flt = DerivativeFilter(gain)
    wvf_out = flt(wvf)
    @test isapprox(wvf_out.signal, gain * vcat(sig[2] - sig[1], diff(sig)...))

    @testset "Test units and precision types" begin
        for T1 in (Float16, Float32, Float64)
            for T2 in (Float16, Float32, Float64)
                sig = rand(T1, 100)
                wvf = RDWaveform(range(0u"ns", step = 16u"ns", length = length(sig)), sig)
                gain = rand(T1) * u"keV"
                flt = DerivativeFilter(gain)
                wvf_out = flt(wvf)
                @test isapprox(wvf_out.signal, gain * vcat(sig[2] - sig[1], diff(sig)...))
            end
        end
    end
end