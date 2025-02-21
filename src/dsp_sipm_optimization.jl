# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_sg_sipm_thresholds_compressed(wvfs::ArrayOfRDWaveforms, config::PropDict)

This function calculates the baseline of the waveforms and the baseline of the waveforms with the sign flipped.
The function is used to calculate the thresholds for the SiPMs.

# Arguments
- `wvfs::ArrayOfRDWaveforms`: Array of RDWaveforms
- `config::PropDict`: Configuration parameters

# Returns
- `Table`: Table with the baseline of the waveforms and the baseline of the waveforms with the sign flipped
"""
function dsp_sg_sipm_thresholds_compressed(wvfs::ArrayOfRDWaveforms, sg_window_length::Unitful.Time{<:Real}, config::PropDict)
    # get dsp meta parameters
    sg_flt_degree = config.sg_flt_degree

    # get waveform data 
    wvfs = decode_data(wvfs)

    # shift waveform by 0 to get Float64 conversation --> ToDO: check if this is necessary
    wvfs = shift_waveform.(wvfs, 0.0)

    # savitzky golay filter: takes derivative of waveform plus smoothing
    sgflt_savitz = SavitzkyGolayFilter(sg_window_length, sg_flt_degree, 1)
    wvfs = sgflt_savitz.(wvfs)

    # project waveforms on the y-axis
    bsl_deriv = vec(flatview(wvfs.signal))

    # remove discharges
    # integrate derivative
    integrator_filter = IntegratorFilter(gain=1)
    wvfs = integrator_filter.(wvfs)

    # project waveforms on the y-axis
    bsl = vec(flatview(wvfs.signal))

    # output Table 
    return TypedTables.Table(
        bsl_deriv = bsl_deriv,
        bsl = bsl,
        bsl_flipped = -bsl
    )
end
export dsp_sg_sipm_thresholds_compressed


"""
    dsp_sg_sipm_optimization_compressed(wvfs::ArrayOfRDWaveforms, dsp_config::PropDict, optimization_config::PropDict)

This function calculates DSP grid to find the optimal thresholds for the SiPMs.

# Arguments
- `wvfs::ArrayOfRDWaveforms`: Array of RDWaveforms
- `dsp_config::PropDict`: Configuration parameters for the DSP
- `optimization_config::PropDict`: Configuration parameters for the optimization

# Returns
- `Table`: Table with the optimal thresholds for the SiPMs
"""
function dsp_sg_sipm_optimization_compressed end
export dsp_sg_sipm_optimization_compressed

function dsp_sg_sipm_optimization_compressed(wvfs::ArrayOfRDWaveforms, dsp_config::PropDict, optimization_config::PropDict)
    # get dsp meta parameters
    min_tot_intersect    = dsp_config.min_tot_intersect
    max_tot_intersect    = dsp_config.max_tot_intersect
    n_σ_threshold        = dsp_config.n_σ_threshold
    sg_flt_degree        = dsp_config.sg_flt_degree
    e_grid_wl            = optimization_config.e_grid_wl
    min_cut_threshold    = optimization_config.threshold.min_cut
    max_cut_threshold    = optimization_config.threshold.max_cut
    
    n_wvfs_threshold = min(length(wvfs), optimization_config.threshold.n_wvfs)  
    
    # get waveform data 
    wvfs = decode_data(wvfs)

    # shift waveform by 0 to get Float64 conversation --> ToDO: check if this is necessary
    wvfs = shift_waveform.(wvfs, 0.0)

    # extract trig max for each waveform based on the grid
    trig_max_grid = Vector{Vector{Float64}}(undef, length(e_grid_wl))
    thresholds_grid = Vector{Float64}(undef, length(e_grid_wl))
    for w in eachindex(e_grid_wl)
        wl = e_grid_wl[w]

        # savitzky golay filter: takes derivative of waveform plus smoothing
        sgflt_savitz = SavitzkyGolayFilter(wl, sg_flt_degree, 1)
        wvfs_sgflt_savitz = sgflt_savitz.(wvfs)

        # get baseline waveforms
        bsl = vec(flatview(wvfs_sgflt_savitz[1:n_wvfs_threshold].signal))
        threshold = thresholdstats(bsl, min_cut_threshold, max_cut_threshold) * n_σ_threshold  
        
        # maximum finder
        intflt = IntersectMaximum(min_tot_intersect, max_tot_intersect)
        inters = intflt.(wvfs_sgflt_savitz, threshold)
        
        trig_max_grid[w] = reduce(vcat, inters.max)
        thresholds_grid[w]  = threshold
    end

    # output Table
    return TypedTables.Table(
        trig_max_grid = VectorOfVectors(trig_max_grid),
        thresholds_grid = thresholds_grid
    )
end

function dsp_sg_sipm_optimization_compressed(n_max_wvfs::Int, wvfs::ArrayOfRDWaveforms, dsp_config::PropDict, optimization_config::PropDict)
    # get dsp meta parameters
    e_grid_wl = optimization_config.e_grid_wl

    # flatten over waveform partitions
    dsp_grid = LegendDataTypes.fast_flatten([let dsp_part = dsp_sg_sipm_optimization_compressed(wf, dsp_config, optimization_config)
        TypedTables.Table(merge(NamedTuple{Tuple(Symbol.(e_grid_wl))}(dsp_part.trig_max_grid...), (thresholds = dsp_part.thresholds_grid, )))
    end for wf in Iterators.partition(decode_data(wvfs), n_max_wvfs)])

    trig_max_grid = Vector{Vector{Float64}}(undef, length(e_grid_wl))
    for (w, wl) in enumerate(e_grid_wl)
        trig_max_grid[w] = getproperty(dsp_grid, Symbol(wl))
    end
    thresholds_grid = [minimum(dsp_grid.thresholds[i:length(e_grid_wl):end]) for i in eachindex(e_grid_wl)]
    return TypedTables.Table(
        trig_max_grid = VectorOfVectors(trig_max_grid),
        thresholds_grid = thresholds_grid
    )
end
