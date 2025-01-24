# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

using Test
using LegendDSP
using RadiationDetectorSignals
using Unitful

@testset "TimeAxisFilter" begin
    @testset "Single waveform" begin
        # create random waveform with random start point
        old_offset = rand()u"ns"
        old_step   = rand()u"ns"
        t = range(old_offset, step = old_step, length = 100)
        signal = rand(100)
        wf = RDWaveform(t, signal)
        @test first(wf.time) == old_offset
        @test step(wf.time)  == old_step
        
        # overwrite step and add offset (using random values to test flexibility)
        new_offset = rand()u"ns"
        new_step   = rand()u"ns"
        flt = TimeAxisFilter(new_step, new_offset)
        wf_new = flt(wf)
        @test first(wf_new.time) == old_offset + new_offset
        @test step(wf_new.time)  == new_step
        
        # explicitly set the start point to zero and use a step of 4ns
        fixed_step = 4.0u"ns"
        flt = TimeAxisFilter(fixed_step, -first(wf.time))
        wf_to0 = flt(wf)
        @test iszero(first(wf_to0.time))
        @test step(wf_to0.time) == fixed_step
    end
    
    @testset "Multiple waveforms" begin
        # create 10 random waveforms with a step of 16ns
        old_offset = 0.0u"ns"
        old_step   = 16.0u"ns"
        t = range(old_offset, step = old_step, length = 100)
        wfs = ArrayOfRDWaveforms(RDWaveform.(Ref(t), rand.(fill(100,10))))
        @test length(wfs) == 10
        @test all(wfs.time .== Ref(t))
        
        new_offset = 100u"ns"
        new_step   =   4u"ns"
        flt = TimeAxisFilter(new_step, new_offset)
        wfs_new = flt.(wfs)
        # the broadcasting result and individual result should be identical
        @test first(wfs_new) == flt(first(wfs))
        # all time axes should be the same
        new_t = flt(first(wfs)).time
        @test all(wfs_new.time .== Ref(new_t))
        @test first(new_t) == old_offset + new_offset
        @test step(new_t) == new_step
    end
end