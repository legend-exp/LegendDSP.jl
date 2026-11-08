using Test
using LegendDSP
using RadiationDetectorDSP
using RadiationDetectorSignals
using Unitful

@testset "Test WeightedSavitzkyGolayFilter" begin
    n = 600
    reltol = 1e-6
    t = range(0u"μs", 20u"μs", 2*n)
    signal = vcat(zeros(n), 10*ones(n))
    wf = RDWaveform(t, signal)

    # define filter parameters and filter
    flt = WeightedSavitzkyGolayFilter(length=1u"μs", degree=3, weightType=2)

    # apply filter to signal
    @test_nowarn flt(wf)

    x = vcat(zeros(10), 10*ones(10))
    flt = WeightedSavitzkyGolayFilter(5, 4, 2)
    res = [
        0.0, 
        0.0, 
        0.0, 
        0.0, 
        0.0, 
        0.0, 
        0.0, 
        0.0, 
        1.4807384272286174e-15, 
        1.592621449357092e-15, 
        9.999999999999998, 
        9.999999999999996, 
        9.999999999999998, 
        9.999999999999998, 
        9.999999999999998, 
        9.999999999999998, 
        9.999999999999998, 
        9.999999999999998, 
        9.999999999999998, 
        10.000000000000002
    ]
    @test isapprox(res, flt(x), rtol=1e-6)   
end

@testset "Test ModifiedSincFilter" begin
    n = 600
    reltol = 1e-6
    t = range(0u"μs", 20u"μs", 2*n)
    signal = vcat(zeros(n), 10*ones(n))
    wf = RDWaveform(t, signal)

    # define filter parameters and filter
    flt = ModifiedSincFilter(d=4, m=1u"μs")

    # apply filter to signal
    @test_nowarn flt(wf)

    x = vcat(zeros(10), 10*ones(10))
    flt = ModifiedSincFilter(2, 3)
    res = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.055151876587386656,
        -0.49321625757844534,
        1.243695939581829,
        8.75630406041817,
        10.493216257578444,
        9.944848123412612,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.000000000000002
    ]
    @test isapprox(res, flt(x); rtol=1e-6)

    # symmetric kernel: the step response is antisymmetric about the step
    y = flt(x)
    @test isapprox(y .+ reverse(y), fill(10.0, length(y)); rtol=1e-6)

    # normalized kernel: a constant passes unchanged, for every degree
    for (d, m) in ((2, 4), (4, 4), (6, 5), (8, 6), (10, 7))
        @test isapprox(ModifiedSincFilter(d, m)(fill(3.0, 30)), fill(3.0, 30); rtol=1e-6)
    end
end

@testset "Test WhittakerHendersonFilter" begin
    n = 600
    reltol = 1e-6
    t = range(0u"μs", 20u"μs", 2*n)
    signal = vcat(zeros(n), 10*ones(n))
    wf = RDWaveform(t, signal)

    # define filter parameters and filter
    flt = WhittakerHendersonFilter(p=3, λ=4)

    # apply filter to signal
    @test_nowarn flt(wf)

    x = vcat(zeros(10), 10*ones(10))
    flt = WhittakerHendersonFilter(p = 3)
    res = [
        -0.02773818585540094,
        0.012639233969518396,
        0.06203483300722407,
        0.09271042540231514,
        0.03409068755810708,
        -0.18030096074202678,
        -0.4731960498109986,
        -0.4498440659762424,
        0.5954235817600229,
        3.199433545249531,
        6.800566454750464,
        9.404576418239975,
        10.449844065976244,
        10.473196049810994,
        10.180300960742013,
        9.96590931244188,
        9.907289574597682,
        9.937965166992784,
        9.987360766030491,
        10.0277381858554
    ]
    @test isapprox(res, flt(x); rtol=1e-6)
end
