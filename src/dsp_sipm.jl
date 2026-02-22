# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    _apply_sipm_filter_dsp(wvfs, flt_config, filter_type::Symbol, wl)

Helper function to apply a single filter configuration and return trigger results.
Returns NamedTuple with trig_max, trig_pos_1, trig_pos_2, trig_pos_diff, and DC versions.
"""
function _apply_sipm_filter_dsp(wvfs_input, flt_config, filter_type::Symbol, wl)
    # Extract filter parameters
    sg_flt_degree    = flt_config.sg_flt_degree
    min_tot_intersect = flt_config.min_tot_intersect
    max_tot_intersect = flt_config.max_tot_intersect
    min_threshold    = flt_config.min_threshold
    max_threshold    = flt_config.max_threshold
    n_σ_threshold    = flt_config.n_σ_threshold
    min_dc_threshold = flt_config.min_dc_threshold
    max_dc_threshold = flt_config.max_dc_threshold
    n_σ_dc_threshold = flt_config.n_σ_dc_threshold
    threshold_method = get(flt_config, :threshold_method, "band_std")
    
    # Create filter based on type
    if filter_type == :sg
        flt = SavitzkyGolayFilter(wl, sg_flt_degree, 1)  # derivative=1
    elseif filter_type == :sgw
        weight_type = haskey(flt_config, :weight_type) ? flt_config.weight_type : 2
        flt = WeightedSavitzkyGolayFilter(length=wl, degree=sg_flt_degree, weightType=weight_type, derivative=1)
    else
        error("Unknown filter type: $filter_type")
    end
    
    # Apply filter
    wvfs_flt = flt.(wvfs_input)
    
    # Maximum finder
    intflt = IntersectMaximum(min_tot_intersect, max_tot_intersect)
    # Threshold calculation: MAD or band_std
    if threshold_method == "mad"
        inters_thres = thresholdstats_mad.(wvfs_flt, min_threshold, max_threshold)
    else
        inters_thres = thresholdstats.(wvfs_flt, min_threshold, max_threshold)
    end
    inters = intflt.(wvfs_flt, n_σ_threshold .* inters_thres)
    
    # For DC detection: integrate and flip
    integrator_filter = IntegratorFilter(gain=1)
    wvfs_int = integrator_filter.(wvfs_flt)
    wvfs_dc = multiply_waveform.(wvfs_int, -1.0)
    if threshold_method == "mad"
        inters_thres_DC = thresholdstats_mad.(wvfs_dc, min_dc_threshold, max_dc_threshold)
    else
        inters_thres_DC = thresholdstats.(wvfs_dc, min_dc_threshold, max_dc_threshold)
    end
    inters_DC = intflt.(wvfs_dc, n_σ_dc_threshold .* inters_thres_DC)
    
    return (
        trig_max = VectorOfVectors(inters.max),
        trig_pos_1 = VectorOfVectors(inters.x),
        trig_pos_2 = VectorOfVectors(inters.x2),
        trig_pos_diff = VectorOfVectors(inters.x_diff),
        trig_max_DC = VectorOfVectors(inters_DC.max),
        trig_pos_1_DC = VectorOfVectors(inters_DC.x),
        trig_pos_2_DC = VectorOfVectors(inters_DC.x2),
        trig_pos_diff_DC = VectorOfVectors(inters_DC.x_diff),
        threshold = inters_thres,
        threshold_DC = inters_thres_DC
    )
end


"""
    _apply_sipm_trap_dsp(wvfs_input, flt_config)

Helper function to apply trap filter directly on raw waveform for SiPM:
1. Trapezoidal filter (directly on raw waveform)
2. Threshold calculation via MAD in [min_threshold, max_threshold] interval
3. IntersectMaximum for trigger finding

Returns NamedTuple with trig_max_trap, trig_pos_1_trap, trig_pos_2_trap, trig_pos_diff_trap, threshold_trap.
"""
function _apply_sipm_trap_dsp(wvfs_input, flt_config)
    # Extract filter parameters
    trap_rt           = flt_config.rt
    trap_ft           = flt_config.ft
    min_tot_intersect = flt_config.min_tot_intersect
    max_tot_intersect = flt_config.max_tot_intersect
    min_threshold     = flt_config.min_threshold
    max_threshold     = flt_config.max_threshold
    n_σ_threshold     = flt_config.n_σ_threshold
    threshold_method  = get(flt_config, :threshold_method, "mad")
    
    # Step 1: Trapezoidal filter directly on raw waveform
    trap_flt = TrapezoidalChargeFilter(trap_rt, trap_ft)
    wvfs_trap = trap_flt.(wvfs_input)
    
    # Step 2: Threshold calculation in [min_threshold, max_threshold] interval
    intflt = IntersectMaximum(min_tot_intersect, max_tot_intersect)
    if threshold_method == "mad"
        inters_thres = thresholdstats_mad.(wvfs_trap, min_threshold, max_threshold)
    else
        inters_thres = thresholdstats.(wvfs_trap, min_threshold, max_threshold)
    end
    
    # Step 3: Trigger finding with n_σ * σ threshold
    inters = intflt.(wvfs_trap, n_σ_threshold .* inters_thres)
    
    return (
        trig_max_trap = VectorOfVectors(inters.max),
        trig_pos_1_trap = VectorOfVectors(inters.x),
        trig_pos_2_trap = VectorOfVectors(inters.x2),
        trig_pos_diff_trap = VectorOfVectors(inters.x_diff),
        threshold_trap = inters_thres
    )
end


"""
    _apply_sipm_trap_sg_dsp(wvfs_input, trap_config, sg_config, sg_wl)

Helper function to apply SG derivative + integration + PZ + trap filter chain for SiPM:
1. Savitzky-Golay filter (derivative) with window length sg_wl
2. Integration of derivative
3. DC tagging on integrated waveform (flip + threshold)
4. Pole-zero correction with configurable τ
5. Trapezoidal filter
6. Threshold calculation via MAD
7. IntersectMaximum for trigger finding

Returns NamedTuple with trig_max_trap_sg, trig_pos_1_trap_sg, DC columns, etc.
"""
function _apply_sipm_trap_sg_dsp(wvfs_input, trap_config, sg_config, sg_wl)
    # Extract trap filter parameters (shared with pure trap)
    trap_rt           = trap_config.rt
    trap_ft           = trap_config.ft
    min_tot_intersect = trap_config.min_tot_intersect
    max_tot_intersect = trap_config.max_tot_intersect
    min_threshold     = trap_config.min_threshold
    max_threshold     = trap_config.max_threshold
    n_σ_threshold     = trap_config.n_σ_threshold
    threshold_method  = get(trap_config, :threshold_method, "mad")
    
    # Extract trap_sg specific parameters
    pz_tau            = get(trap_config, :pz_tau, 3.0u"µs")
    n_σ_dc_threshold  = get(trap_config, :n_σ_dc_threshold, 5.0)
    min_dc_threshold  = get(trap_config, :min_dc_threshold, -3.0)
    max_dc_threshold  = get(trap_config, :max_dc_threshold, 3.0)
    
    # Extract SG filter parameters
    sg_flt_degree = sg_config.sg_flt_degree
    
    # Step 1: Savitzky-Golay filter (derivative)
    sg_flt = SavitzkyGolayFilter(sg_wl, sg_flt_degree, 1)  # derivative=1
    wvfs_sg = sg_flt.(wvfs_input)
    
    # Step 2: Integration of derivative
    integrator_filter = IntegratorFilter(gain=1)
    wvfs_int = integrator_filter.(wvfs_sg)
    
    # Step 3: DC tagging on integrated waveform (flip and threshold)
    wvfs_dc = multiply_waveform.(wvfs_int, -1.0)
    intflt_dc = IntersectMaximum(min_tot_intersect, max_tot_intersect)
    if threshold_method == "mad"
        inters_thres_dc = thresholdstats_mad.(wvfs_dc, min_dc_threshold, max_dc_threshold)
    else
        inters_thres_dc = thresholdstats.(wvfs_dc, min_dc_threshold, max_dc_threshold)
    end
    inters_dc = intflt_dc.(wvfs_dc, n_σ_dc_threshold .* inters_thres_dc)
    
    # Step 4: Pole-zero correction (InvCRFilter compensates preamp exponential decay)
    pz_flt = InvCRFilter(pz_tau)
    wvfs_pz = pz_flt.(wvfs_int)
    
    # Step 5: Trapezoidal filter
    trap_flt = TrapezoidalChargeFilter(trap_rt, trap_ft)
    wvfs_trap = trap_flt.(wvfs_pz)
    
    # Step 6: Threshold calculation in [min_threshold, max_threshold] interval
    intflt = IntersectMaximum(min_tot_intersect, max_tot_intersect)
    if threshold_method == "mad"
        inters_thres = thresholdstats_mad.(wvfs_trap, min_threshold, max_threshold)
    else
        inters_thres = thresholdstats.(wvfs_trap, min_threshold, max_threshold)
    end
    
    # Step 7: Trigger finding with n_σ * σ threshold
    inters = intflt.(wvfs_trap, n_σ_threshold .* inters_thres)
    
    return (
        trig_max_trap_sg = VectorOfVectors(inters.max),
        trig_pos_1_trap_sg = VectorOfVectors(inters.x),
        trig_pos_2_trap_sg = VectorOfVectors(inters.x2),
        trig_pos_diff_trap_sg = VectorOfVectors(inters.x_diff),
        threshold_trap_sg = inters_thres,
        # DC tagging columns
        trig_max_DC_trap_sg = VectorOfVectors(inters_dc.max),
        trig_pos_1_DC_trap_sg = VectorOfVectors(inters_dc.x),
        trig_pos_2_DC_trap_sg = VectorOfVectors(inters_dc.x2),
        trig_pos_diff_DC_trap_sg = VectorOfVectors(inters_dc.x_diff),
        threshold_DC_trap_sg = inters_thres_dc
    )
end

"""
    _build_filter_output_names(filter_type::Symbol, wl_ns::Int)

Build column names for a filter configuration, e.g., trig_max_sg_100.
"""
function _build_filter_output_names(filter_type::Symbol, wl_ns::Int)
    suffix = "$(filter_type)_$(wl_ns)"
    return (
        trig_max = Symbol("trig_max_$suffix"),
        trig_pos_1 = Symbol("trig_pos_1_$suffix"),
        trig_pos_2 = Symbol("trig_pos_2_$suffix"),
        trig_pos_diff = Symbol("trig_pos_diff_$suffix"),
        trig_max_DC = Symbol("trig_max_DC_$suffix"),
        trig_pos_1_DC = Symbol("trig_pos_1_DC_$suffix"),
        trig_pos_2_DC = Symbol("trig_pos_2_DC_$suffix"),
        trig_pos_diff_DC = Symbol("trig_pos_diff_DC_$suffix"),
        threshold = Symbol("threshold_$suffix"),
        threshold_DC = Symbol("threshold_DC_$suffix")
    )
end

"""
    dsp_sipm(data::Q, config::PropDict; sg_wl_trap=nothing)

DSP routine for SiPM data with multiple filter configurations (sg/sgw with various window lengths).

# Input data
- `waveform`: waveform data
- `baseline`: baseline data
- `timestamp`: timestamp data  
- `eventnumber`: event number data
- `daqenergy`: energy data

# Output data
Base columns: blfc, timestamp, eventID_fadc, e_fc, t_max, t_min, t_max_lar, t_min_lar, e_max, e_min, e_max_lar, e_min_lar

Per filter/window (e.g., sg_100, sg_150, sg_200, sgw_100, sgw_150, sgw_200):
- trig_max_{filter}_{wl}: trigger maxima
- trig_pos_1_{filter}_{wl}: trigger positions (up-crossing)
- trig_pos_2_{filter}_{wl}: trigger positions (down-crossing)
- trig_pos_diff_{filter}_{wl}: time-over-threshold (trig_pos_2 - trig_pos_1)
- trig_max_DC_{filter}_{wl}: DC trigger maxima (on flipped integrated waveform)
- trig_pos_1_DC_{filter}_{wl}: DC trigger positions (up-crossing)
- trig_pos_2_DC_{filter}_{wl}: DC trigger positions (down-crossing)
- trig_pos_diff_DC_{filter}_{wl}: DC time-over-threshold

If trap filter config exists:
- trig_max_trap: trigger maxima from direct trap filter on raw waveform
- trig_pos_1_trap: trigger positions (up-crossing) from trap filter
- trig_pos_2_trap: trigger positions (down-crossing) from trap filter
- trig_pos_diff_trap: time-over-threshold from trap filter
- threshold_trap: threshold used for trap filter

If trap_sg filter is enabled (sg_wl_trap parameter or trap.sg_wl in config):
- trig_max_trap_sg: trigger maxima from SG+Integration+PZ+Trap chain
- trig_pos_1_trap_sg: trigger positions (up-crossing) from trap_sg
- trig_pos_2_trap_sg: trigger positions (down-crossing) from trap_sg
- trig_pos_diff_trap_sg: time-over-threshold from trap_sg
- threshold_trap_sg: threshold used for trap_sg
- trig_max_DC_trap_sg: DC trigger maxima (on flipped integrated waveform before PZ)
- trig_pos_1_DC_trap_sg, trig_pos_2_DC_trap_sg, trig_pos_diff_DC_trap_sg, threshold_DC_trap_sg

# Config parameters for trap_sg (in filters.trap section):
- sg_wl: SG window length for trap_sg (default from sg_wl_trap parameter)
- pz_tau: Pole-zero time constant (default 3.0µs)
- n_σ_dc_threshold: DC threshold multiplier (default 5.0)
- min_dc_threshold, max_dc_threshold: DC baseline window (default -3.0, 3.0)
"""
function dsp_sipm(data::Q, config::PropDict; sg_wl_trap=nothing) where {Q <: Table}
    # Get common parameters
    t0_hpge_window = config.t0_hpge_window
    
    # Get waveform data 
    wvfs = data.waveform
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # Convert to Float64
    wvfs = shift_waveform.(wvfs, 0.0)

    # Get wvf extrema
    estats = extremestats.(wvfs)
    
    # Get wvf extrema for truncated waveform (LAr window)
    uflt_trunc = TruncateFilter(first(t0_hpge_window)..last(t0_hpge_window))
    estats_trunc = extremestats.(uflt_trunc.(wvfs))

    # Build base output
    output = Dict{Symbol, Any}(
        :blfc => blfc,
        :timestamp => ts,
        :eventID_fadc => evID,
        :e_fc => efc,
        :t_max => uconvert.(u"µs", estats.tmax),
        :t_min => uconvert.(u"µs", estats.tmin),
        :t_max_lar => uconvert.(u"µs", estats_trunc.tmax),
        :t_min_lar => uconvert.(u"µs", estats_trunc.tmin),
        :e_max => estats.max,
        :e_min => estats.min,
        :e_max_lar => estats_trunc.max,
        :e_min_lar => estats_trunc.min
    )
    
    # Process each filter type
    filters_config = config.filters
    for (filter_type, flt_config) in pairs(filters_config)
        filter_sym = Symbol(filter_type)
        
        # Skip trap filter here - handled separately below
        if filter_sym == :trap
            continue
        end
        
        window_lengths = flt_config.window_lengths
        
        for wl in window_lengths
            # Get window length in ns for naming
            wl_ns = round(Int, ustrip(u"ns", wl))
            
            # Apply filter DSP
            result = _apply_sipm_filter_dsp(wvfs, flt_config, filter_sym, wl)
            
            # Build column names
            names = _build_filter_output_names(filter_sym, wl_ns)
            
            # Add results to output
            output[names.trig_max] = result.trig_max
            output[names.trig_pos_1] = result.trig_pos_1
            output[names.trig_pos_2] = result.trig_pos_2
            output[names.trig_pos_diff] = result.trig_pos_diff
            output[names.trig_max_DC] = result.trig_max_DC
            output[names.trig_pos_1_DC] = result.trig_pos_1_DC
            output[names.trig_pos_2_DC] = result.trig_pos_2_DC
            output[names.trig_pos_diff_DC] = result.trig_pos_diff_DC
            output[names.threshold] = result.threshold
            output[names.threshold_DC] = result.threshold_DC
        end
    end
    
    # Process trap filters if trap config exists
    if haskey(filters_config, :trap)
        trap_config = filters_config.trap
        
        # Variant 1: Direct trap filter on raw waveform (always applied if trap config exists)
        trap_result = _apply_sipm_trap_dsp(wvfs, trap_config)
        output[:trig_max_trap] = trap_result.trig_max_trap
        output[:trig_pos_1_trap] = trap_result.trig_pos_1_trap
        output[:trig_pos_2_trap] = trap_result.trig_pos_2_trap
        output[:trig_pos_diff_trap] = trap_result.trig_pos_diff_trap
        output[:threshold_trap] = trap_result.threshold_trap
        
        # Variant 2: SG derivative + integration + PZ + trap filter (with DC tagging)
        # Use sg_wl_trap parameter if provided, otherwise fall back to trap.sg_wl from config
        _sg_wl_trap = !isnothing(sg_wl_trap) ? sg_wl_trap : get(trap_config, :sg_wl, nothing)
        if !isnothing(_sg_wl_trap) && haskey(filters_config, :sg)
            sg_config = filters_config.sg
            trap_sg_result = _apply_sipm_trap_sg_dsp(wvfs, trap_config, sg_config, _sg_wl_trap)
            output[:trig_max_trap_sg] = trap_sg_result.trig_max_trap_sg
            output[:trig_pos_1_trap_sg] = trap_sg_result.trig_pos_1_trap_sg
            output[:trig_pos_2_trap_sg] = trap_sg_result.trig_pos_2_trap_sg
            output[:trig_pos_diff_trap_sg] = trap_sg_result.trig_pos_diff_trap_sg
            output[:threshold_trap_sg] = trap_sg_result.threshold_trap_sg
            # DC tagging columns for trap_sg
            output[:trig_max_DC_trap_sg] = trap_sg_result.trig_max_DC_trap_sg
            output[:trig_pos_1_DC_trap_sg] = trap_sg_result.trig_pos_1_DC_trap_sg
            output[:trig_pos_2_DC_trap_sg] = trap_sg_result.trig_pos_2_DC_trap_sg
            output[:trig_pos_diff_DC_trap_sg] = trap_sg_result.trig_pos_diff_DC_trap_sg
            output[:threshold_DC_trap_sg] = trap_sg_result.threshold_DC_trap_sg
        end
    end
    
    return TypedTables.Table(NamedTuple(output))
end
export dsp_sipm


"""
    dsp_sipm_compressed(data::Q, config::PropDict; sg_wl_trap=nothing)

DSP routine for compressed SiPM data with multiple filter configurations.
Same as dsp_sipm but uses compressed waveform data (waveform_bit_drop).

See `dsp_sipm` for full documentation of output columns and config parameters.
"""
function dsp_sipm_compressed(data::Q, config::PropDict; sg_wl_trap=nothing) where {Q <: Table}
    # Get common parameters
    t0_hpge_window = config.t0_hpge_window
    
    # Get waveform data (compressed)
    wvfs = decode_data(data.waveform_bit_drop)
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # Convert to Float64
    wvfs = shift_waveform.(wvfs, 0.0)

    # Get wvf extrema
    estats = extremestats.(wvfs)
    
    # Get wvf extrema for truncated waveform (LAr window)
    uflt_trunc = TruncateFilter(first(t0_hpge_window)..last(t0_hpge_window))
    estats_trunc = extremestats.(uflt_trunc.(wvfs))

    # Build base output
    output = Dict{Symbol, Any}(
        :blfc => blfc,
        :timestamp => ts,
        :eventID_fadc => evID,
        :e_fc => efc,
        :t_max => uconvert.(u"µs", estats.tmax),
        :t_min => uconvert.(u"µs", estats.tmin),
        :t_max_lar => uconvert.(u"µs", estats_trunc.tmax),
        :t_min_lar => uconvert.(u"µs", estats_trunc.tmin),
        :e_max => estats.max,
        :e_min => estats.min,
        :e_max_lar => estats_trunc.max,
        :e_min_lar => estats_trunc.min
    )
    
    # Process each filter type
    filters_config = config.filters
    for (filter_type, flt_config) in pairs(filters_config)
        filter_sym = Symbol(filter_type)
        
        # Skip trap filter here - handled separately below
        if filter_sym == :trap
            continue
        end
        
        window_lengths = flt_config.window_lengths
        
        for wl in window_lengths
            # Get window length in ns for naming
            wl_ns = round(Int, ustrip(u"ns", wl))
            
            # Apply filter DSP
            result = _apply_sipm_filter_dsp(wvfs, flt_config, filter_sym, wl)
            
            # Build column names
            names = _build_filter_output_names(filter_sym, wl_ns)
            
            # Add results to output
            output[names.trig_max] = result.trig_max
            output[names.trig_pos_1] = result.trig_pos_1
            output[names.trig_pos_2] = result.trig_pos_2
            output[names.trig_pos_diff] = result.trig_pos_diff
            output[names.trig_max_DC] = result.trig_max_DC
            output[names.trig_pos_1_DC] = result.trig_pos_1_DC
            output[names.trig_pos_2_DC] = result.trig_pos_2_DC
            output[names.trig_pos_diff_DC] = result.trig_pos_diff_DC
            output[names.threshold] = result.threshold
            output[names.threshold_DC] = result.threshold_DC
        end
    end
    
    # Process trap filters if trap config exists
    if haskey(filters_config, :trap)
        trap_config = filters_config.trap
        
        # Variant 1: Direct trap filter on raw waveform (always applied if trap config exists)
        trap_result = _apply_sipm_trap_dsp(wvfs, trap_config)
        output[:trig_max_trap] = trap_result.trig_max_trap
        output[:trig_pos_1_trap] = trap_result.trig_pos_1_trap
        output[:trig_pos_2_trap] = trap_result.trig_pos_2_trap
        output[:trig_pos_diff_trap] = trap_result.trig_pos_diff_trap
        output[:threshold_trap] = trap_result.threshold_trap
        
        # Variant 2: SG derivative + integration + PZ + trap filter (with DC tagging)
        # Use sg_wl_trap parameter if provided, otherwise fall back to trap.sg_wl from config
        _sg_wl_trap = !isnothing(sg_wl_trap) ? sg_wl_trap : get(trap_config, :sg_wl, nothing)
        if !isnothing(_sg_wl_trap) && haskey(filters_config, :sg)
            sg_config = filters_config.sg
            trap_sg_result = _apply_sipm_trap_sg_dsp(wvfs, trap_config, sg_config, _sg_wl_trap)
            output[:trig_max_trap_sg] = trap_sg_result.trig_max_trap_sg
            output[:trig_pos_1_trap_sg] = trap_sg_result.trig_pos_1_trap_sg
            output[:trig_pos_2_trap_sg] = trap_sg_result.trig_pos_2_trap_sg
            output[:trig_pos_diff_trap_sg] = trap_sg_result.trig_pos_diff_trap_sg
            output[:threshold_trap_sg] = trap_sg_result.threshold_trap_sg
            # DC tagging columns for trap_sg
            output[:trig_max_DC_trap_sg] = trap_sg_result.trig_max_DC_trap_sg
            output[:trig_pos_1_DC_trap_sg] = trap_sg_result.trig_pos_1_DC_trap_sg
            output[:trig_pos_2_DC_trap_sg] = trap_sg_result.trig_pos_2_DC_trap_sg
            output[:trig_pos_diff_DC_trap_sg] = trap_sg_result.trig_pos_diff_DC_trap_sg
            output[:threshold_DC_trap_sg] = trap_sg_result.threshold_DC_trap_sg
        end
    end
    
    return TypedTables.Table(NamedTuple(output))
end
export dsp_sipm_compressed
