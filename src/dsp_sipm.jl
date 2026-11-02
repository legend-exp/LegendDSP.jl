# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    dsp_sipm(data::Q, config::PropDict, pars::PropDict)

DSP routine for SiPM data. `pars` carries the shaping lengths of the optimisation (`sg.wl`, `mwa.wl`) and, when the
template step has run for the channel, its inverse `kernel` (`taps`, `dt`). The MWA chain (`config.filters.mwa`:
`wl_trigger`, `passes`, `n_σ_threshold`, `min_threshold`, `max_threshold`, `solve_overlap`, `merge_distance`) sends
the waveform through the kernel (the derivative without it), shapes it with `passes` moving averages of the trigger
length and of the energy length `mwa.wl`, both scaled so that a unit stick (one photon) has peak 1, and finds the
sites. With `solve_overlap` the sites are the extrema of the trigger trace and their amplitudes are solved on the
energy trace (matrix method); without it, sites and amplitudes are the extrema of the energy trace. Sites within
`2 passes` times the longer shaping length of either trace end are dropped.

# Input data
The input data is a table with the following columns:
- `waveform`: waveform data
- `baseline`: baseline data
- `timestamp`: timestamp data
- `eventnumber`: event number data
- `daqenergy`: energy data

# Output data
The output data is a table with the following columns:
- `blfc`: baseline from FADC
- `timestamp`: timestamp
- `eventID_fadc`: event number from FADC
- `e_fc`: energy from FADC
- `t_max`, `t_min`: time of waveform maximum/minimum
- `t_max_lar`, `t_min_lar`: time of waveform maximum/minimum in LAr window
- `e_max`, `e_min`: waveform maximum/minimum amplitude
- `e_max_lar`, `e_min_lar`: waveform maximum/minimum amplitude in LAr window
- `blmean`, `blsigma`, `blslope`, `bloffset`: baseline statistics
- `wfmean`, `wfsigma`, `wfslope`, `wfoffset`: waveform statistics
## SG filter columns
- `threshold_sg`: SG trigger threshold
- `threshold_DC_sg`: SG discharge trigger threshold
- `trig_pos_sg`: SG trigger positions
- `trig_max_sg`: SG trigger maxima
- `trig_pos_DC_sg`: SG discharge trigger positions
- `trig_max_DC_sg`: SG discharge trigger maxima
## Trap filter columns
- `threshold_trap`: trap trigger threshold
- `threshold_DC_trap`: trap discharge trigger threshold
- `trig_pos_trap`: trap trigger positions
- `trig_pos_high_trap`: trap trigger high positions
- `trig_pos_tot_trap`: trap trigger time-over-threshold positions
- `trig_max_trap`: trap trigger maxima
- `trig_pos_DC_trap`: trap discharge trigger positions
- `trig_pos_high_DC_trap`: trap discharge trigger high positions
- `trig_pos_tot_DC_trap`: trap discharge trigger time-over-threshold positions
- `trig_max_DC_trap`: trap discharge trigger maxima
## MWA filter columns
- `threshold_mwa`: MAD-σ of the site trace (the trigger trace with the matrix method, the energy trace without)
- `trig_pos_mwa`: site positions, every local maximum above `n_σ_threshold` σ and minimum below its negative
- `trig_pos_tot_mwa`: time over threshold of the excursion a site lies in
- `trig_max_mwa`: site amplitudes on the energy shaping, signed; solved for the overlap with the matrix method
"""
function dsp_sipm(data::Q, config::PropDict, pars::PropDict) where {Q <: Table}
    # common config
    sg_flt_degree         = config.sg_flt_degree
    t0_hpge_window        = config.t0_hpge_window

    # get SG config parameters
    sg_config = config.filters.sg
    min_tot_intersect     = sg_config.min_tot_intersect
    max_tot_intersect     = sg_config.max_tot_intersect
    min_threshold         = sg_config.min_threshold
    max_threshold         = sg_config.max_threshold
    n_σ_threshold         = sg_config.n_σ_threshold
    min_dc_threshold      = sg_config.min_dc_threshold
    max_dc_threshold      = sg_config.max_dc_threshold
    n_σ_dc_threshold      = sg_config.n_σ_dc_threshold

    # get trap config parameters
    trap_config = config.filters.trap
    trap_rt                = trap_config.rt
    trap_ft                = trap_config.ft
    trap_pz_tau            = trap_config.pz_tau
    trap_min_tot_intersect = trap_config.min_tot_intersect
    trap_max_tot_intersect = trap_config.max_tot_intersect
    trap_min_threshold     = trap_config.min_threshold
    trap_max_threshold     = trap_config.max_threshold
    trap_n_σ_threshold     = trap_config.n_σ_threshold
    trap_min_dc_threshold  = trap_config.min_dc_threshold
    trap_max_dc_threshold  = trap_config.max_dc_threshold
    trap_n_σ_dc_threshold  = trap_config.n_σ_dc_threshold

    # get MWA config parameters
    mwa_config = config.filters.mwa
    mwa_wl_trigger         = mwa_config.wl_trigger
    mwa_passes             = mwa_config.passes
    mwa_n_σ_threshold      = mwa_config.n_σ_threshold
    mwa_min_threshold      = mwa_config.min_threshold
    mwa_max_threshold      = mwa_config.max_threshold
    mwa_solve_overlap      = mwa_config.solve_overlap
    mwa_merge_distance     = mwa_config.merge_distance
    iseven(mwa_passes) || throw(ArgumentError("The MWA chain needs an even number of passes to stay centred, got $mwa_passes"))

    # unpack optimization parameters and, if the template step has run, the inverse kernel
    sg_window_length = pars.sg.wl
    mwa_wl_energy    = pars.mwa.wl
    has_kernel       = haskey(pars, :kernel)

    # get waveform data
    data = decode_data(data)
    wvfs = data.waveform
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # shift waveform by 0 to get Float64 conversion
    wvfs = shift_waveform.(wvfs, 0.0)

    # get wvf maximum and minimum with timepoints
    estats = extremestats.(wvfs)

    # get wvf maximum and minimum with timepoints for truncated waveform
    uflt_trunc = TruncateFilter(first(t0_hpge_window)..last(t0_hpge_window))
    estats_trunc = extremestats.(uflt_trunc.(wvfs))

    # === SG pipeline ===
    # savitzky golay filter: takes derivative of waveform plus smoothing
    sgflt_savitz = SavitzkyGolayFilter(sg_window_length, sg_flt_degree, 1)
    wvfs_sg = sgflt_savitz.(wvfs)

    # trigger finding on SG-filtered waveforms
    intflt_sg = IntersectMaximum(min_tot_intersect, max_tot_intersect)
    inters_thres = thresholdstats_mad.(wvfs_sg, min_threshold, max_threshold)
    inters = intflt_sg.(wvfs_sg, n_σ_threshold .* inters_thres)

    # integrate derivative
    integrator_filter = IntegratorFilter(gain=1)
    wvfs_int = integrator_filter.(wvfs_sg)

    # get blstats on integrated waveforms
    time_min = minimum(wvfs_int[1].time)
    time_max = last(wvfs_int[1].time)
    # baseline window ends at the first trigger (end of waveform if none), but spans at least 3 samples
    t_bl_min = time_min + 3*step(wvfs_int[1].time)
    t_bl_stop = map(x -> oftype(time_max, max(minimum(x; init = time_max), t_bl_min)), inters.x)
    bl_stats = signalstats.(wvfs_int, time_min, t_bl_stop)
    sigstats = signalstats.(wvfs_int, time_min, time_max)
    
    # flip around x-axis for discharge detection
    wvfs_int_flip = multiply_waveform.(wvfs_int, -1.0)
    inters_thres_DC = thresholdstats_mad.(wvfs_int_flip, min_dc_threshold, max_dc_threshold)
    inters_DC = intflt_sg.(wvfs_int_flip, n_σ_dc_threshold .* inters_thres_DC)

    # === Trap pipeline ===
    # PZ correction on integrated waveforms
    pz_filter = InvCRFilter(trap_pz_tau)
    wvfs_pz = pz_filter.(wvfs_int)

    # trapezoidal charge filter
    trap_filter = TrapezoidalChargeFilter(trap_rt, trap_ft)
    wvfs_trap = trap_filter.(wvfs_pz)

    # trigger finding on trap-filtered waveforms
    intflt_trap = IntersectMaximum(trap_min_tot_intersect, trap_max_tot_intersect)
    inters_thres_trap = thresholdstats_mad.(wvfs_trap, trap_min_threshold, trap_max_threshold)
    inters_trap = intflt_trap.(wvfs_trap, trap_n_σ_threshold .* inters_thres_trap)

    # discharge detection on flipped integrated waveforms 
    inters_thres_DC_trap = thresholdstats_mad.(wvfs_int_flip, trap_min_dc_threshold, trap_max_dc_threshold)
    inters_DC_trap = intflt_sg.(wvfs_int_flip, trap_n_σ_dc_threshold .* inters_thres_DC_trap)

    # === MWA pipeline ===
    # photon trace: the baseline-subtracted waveform through the inverse kernel of the template step; a channel without
    # template pars takes the derivative (kernel [1, -1]) instead
    t_axis = first(wvfs).time
    dt = step(t_axis)
    !has_kernel || pars.kernel.dt ≈ dt || throw(ArgumentError("The kernel taps are sampled at $(pars.kernel.dt), the waveforms at $dt"))
    photon_filter = CenteredConvolutionFilter(has_kernel ? Float64.(pars.kernel.taps) : [0.0, 1.0, -1.0])
    wvfs_bl = shift_waveform.(wvfs, [-median(wf.signal) for wf in wvfs])
    wvfs_photons = photon_filter.(wvfs_bl)

    # short (trigger) and long (energy) shaping, each scaled so that a unit stick (one photon after the kernel) has peak 1
    mwa_trigger_filter = MovingWindowMultiFilter(mwa_wl_trigger; passes = mwa_passes)
    mwa_energy_filter  = MovingWindowMultiFilter(mwa_wl_energy; passes = mwa_passes)
    stick = RDWaveform(t_axis, [i == length(t_axis) ÷ 2 ? 1.0 : 0.0 for i in eachindex(t_axis)])
    response_trigger = mwa_trigger_filter(stick).signal
    response_energy  = mwa_energy_filter(stick).signal
    wvfs_trigger = multiply_waveform.(mwa_trigger_filter.(wvfs_photons), 1 / maximum(response_trigger))
    wvfs_energy  = multiply_waveform.(mwa_energy_filter.(wvfs_photons), 1 / maximum(response_energy))
    wvfs_trigger = shift_waveform.(wvfs_trigger, [-median(wf.signal) for wf in wvfs_trigger])
    wvfs_energy  = shift_waveform.(wvfs_energy, [-median(wf.signal) for wf in wvfs_energy])

    # sites: the extrema beyond n σ of the trigger trace (matrix method) or of the energy trace
    wvfs_sites = mwa_solve_overlap ? wvfs_trigger : wvfs_energy
    mwa_thresholds = thresholdstats_mad.(wvfs_sites, mwa_min_threshold, mwa_max_threshold)
    sites = extrema_over_threshold.(wvfs_sites, mwa_n_σ_threshold .* mwa_thresholds)

    # matrix method: the amplitudes of all sites solved on the long shaping (otherwise the extrema of the energy trace)
    if mwa_solve_overlap
        solver = SiteSolver(response_energy ./ maximum(response_energy), argmax(response_energy), mwa_merge_distance * uconvert(NoUnits, mwa_wl_energy / dt))
        sites = map(zip(wvfs_energy, sites)) do (wf, s)
            sol = solver(wf.signal, uconvert.(NoUnits, (s.x .- first(t_axis)) ./ dt) .+ 1, uconvert.(NoUnits, s.x_tot ./ dt), Int.(sign.(s.max)))
            (x = convert.(eltype(t_axis), first(t_axis) .+ (sol.pos .- 1) .* dt), x_tot = convert.(eltype(t_axis), sol.tot .* dt), max = sol.A)
        end
    end

    # without the shaping transients at both ends of the trace
    mwa_margin = 2 * mwa_passes * max(mwa_wl_trigger, mwa_wl_energy)
    sites = map(sites) do s
        k = findall(t -> first(t_axis) + mwa_margin <= t <= last(t_axis) - mwa_margin, s.x)
        (x = s.x[k], x_tot = s.x_tot[k], max = s.max[k])
    end

    # output Table 
    return TypedTables.Table(
        blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
        t_max = uconvert.(u"µs", estats.tmax), t_min = uconvert.(u"µs", estats.tmin), t_max_lar = uconvert.(u"µs", estats_trunc.tmax), t_min_lar = uconvert.(u"µs", estats_trunc.tmin),
        e_max = estats.max, e_min = estats.min, e_max_lar = estats_trunc.max, e_min_lar = estats_trunc.min,
        blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
        wfmean = sigstats.mean, wfsigma = sigstats.sigma, wfslope = sigstats.slope, wfoffset = sigstats.offset,
        # SG triggers
        threshold_sg = inters_thres, threshold_DC_sg = inters_thres_DC,
        trig_pos_sg = VectorOfVectors(inters.x), trig_max_sg = VectorOfVectors(inters.max),
        trig_pos_DC_sg = VectorOfVectors(inters_DC.x), trig_max_DC_sg = VectorOfVectors(inters_DC.max),
        # Trap triggers
        threshold_trap = inters_thres_trap, threshold_DC_trap = inters_thres_DC_trap,
        trig_pos_trap = VectorOfVectors(inters_trap.x), trig_pos_high_trap = VectorOfVectors(inters_trap.x_high),
        trig_pos_tot_trap = VectorOfVectors(inters_trap.x_tot), trig_max_trap = VectorOfVectors(inters_trap.max),
        trig_pos_DC_trap = VectorOfVectors(inters_DC_trap.x), trig_pos_high_DC_trap = VectorOfVectors(inters_DC_trap.x_high),
        trig_pos_tot_DC_trap = VectorOfVectors(inters_DC_trap.x_tot), trig_max_DC_trap = VectorOfVectors(inters_DC_trap.max),
        # MWA triggers
        threshold_mwa = mwa_thresholds, trig_pos_mwa = VectorOfVectors([s.x for s in sites]),
        trig_pos_tot_mwa = VectorOfVectors([s.x_tot for s in sites]), trig_max_mwa = VectorOfVectors([s.max for s in sites])
    )
end
export dsp_sipm



"""
    dsp_sipm_compressed(data::Q, config::PropDict, pars::PropDict)

DSP routine for SiPM data, see `dsp_sipm`; reads `waveform_bit_drop`.

# Input data
The input data is a table with the following columns:
- `waveform_bit_drop`: bit-drop compressed waveform data
- `baseline`: baseline data
- `timestamp`: timestamp data
- `eventnumber`: event number data
- `daqenergy`: energy data

# Output data
The output data is a table with the following columns:
- `blfc`: baseline from FADC
- `timestamp`: timestamp
- `eventID_fadc`: event number from FADC
- `e_fc`: energy from FADC
- `t_max`, `t_min`: time of waveform maximum/minimum
- `t_max_lar`, `t_min_lar`: time of waveform maximum/minimum in LAr window
- `e_max`, `e_min`: waveform maximum/minimum amplitude
- `e_max_lar`, `e_min_lar`: waveform maximum/minimum amplitude in LAr window
- `blmean`, `blsigma`, `blslope`, `bloffset`: baseline statistics
- `wfmean`, `wfsigma`, `wfslope`, `wfoffset`: waveform statistics
## SG filter columns
- `threshold_sg`: SG trigger threshold
- `threshold_DC_sg`: SG discharge trigger threshold
- `trig_pos_sg`: SG trigger positions
- `trig_max_sg`: SG trigger maxima
- `trig_pos_DC_sg`: SG discharge trigger positions
- `trig_max_DC_sg`: SG discharge trigger maxima
## Trap filter columns
- `threshold_trap`: trap trigger threshold
- `threshold_DC_trap`: trap discharge trigger threshold
- `trig_pos_trap`: trap trigger positions
- `trig_pos_high_trap`: trap trigger high positions
- `trig_pos_tot_trap`: trap trigger time-over-threshold positions
- `trig_max_trap`: trap trigger maxima
- `trig_pos_DC_trap`: trap discharge trigger positions
- `trig_pos_high_DC_trap`: trap discharge trigger high positions
- `trig_pos_tot_DC_trap`: trap discharge trigger time-over-threshold positions
- `trig_max_DC_trap`: trap discharge trigger maxima
## MWA filter columns
- `threshold_mwa`: MAD-σ of the site trace (the trigger trace with the matrix method, the energy trace without)
- `trig_pos_mwa`: site positions, every local maximum above `n_σ_threshold` σ and minimum below its negative
- `trig_pos_tot_mwa`: time over threshold of the excursion a site lies in
- `trig_max_mwa`: site amplitudes on the energy shaping, signed; solved for the overlap with the matrix method
"""
function dsp_sipm_compressed(data::Q, config::PropDict, pars::PropDict) where {Q <: Table}
    # common config
    sg_flt_degree         = config.sg_flt_degree
    t0_hpge_window        = config.t0_hpge_window

    # get SG config parameters
    sg_config = config.filters.sg
    min_tot_intersect     = sg_config.min_tot_intersect
    max_tot_intersect     = sg_config.max_tot_intersect
    min_threshold         = sg_config.min_threshold
    max_threshold         = sg_config.max_threshold
    n_σ_threshold         = sg_config.n_σ_threshold
    min_dc_threshold      = sg_config.min_dc_threshold
    max_dc_threshold      = sg_config.max_dc_threshold
    n_σ_dc_threshold      = sg_config.n_σ_dc_threshold

    # get trap config parameters
    trap_config = config.filters.trap
    trap_rt                = trap_config.rt
    trap_ft                = trap_config.ft
    trap_pz_tau            = trap_config.pz_tau
    trap_min_tot_intersect = trap_config.min_tot_intersect
    trap_max_tot_intersect = trap_config.max_tot_intersect
    trap_min_threshold     = trap_config.min_threshold
    trap_max_threshold     = trap_config.max_threshold
    trap_n_σ_threshold     = trap_config.n_σ_threshold
    trap_min_dc_threshold  = trap_config.min_dc_threshold
    trap_max_dc_threshold  = trap_config.max_dc_threshold
    trap_n_σ_dc_threshold  = trap_config.n_σ_dc_threshold

    # get MWA config parameters
    mwa_config = config.filters.mwa
    mwa_wl_trigger         = mwa_config.wl_trigger
    mwa_passes             = mwa_config.passes
    mwa_n_σ_threshold      = mwa_config.n_σ_threshold
    mwa_min_threshold      = mwa_config.min_threshold
    mwa_max_threshold      = mwa_config.max_threshold
    mwa_solve_overlap      = mwa_config.solve_overlap
    mwa_merge_distance     = mwa_config.merge_distance
    iseven(mwa_passes) || throw(ArgumentError("The MWA chain needs an even number of passes to stay centred, got $mwa_passes"))

    # unpack optimization parameters and, if the template step has run, the inverse kernel
    sg_window_length = pars.sg.wl
    mwa_wl_energy    = pars.mwa.wl
    has_kernel       = haskey(pars, :kernel)

    # get waveform data
    data = decode_data(data)
    wvfs = data.waveform_bit_drop
    blfc = data.baseline
    ts   = data.timestamp
    evID = data.eventnumber
    efc  = data.daqenergy

    # shift waveform by 0 to get Float64 conversion
    wvfs = shift_waveform.(wvfs, 0.0)

    # get wvf maximum and minimum with timepoints
    estats = extremestats.(wvfs)

    # get wvf maximum and minimum with timepoints for truncated waveform
    uflt_trunc = TruncateFilter(first(t0_hpge_window)..last(t0_hpge_window))
    estats_trunc = extremestats.(uflt_trunc.(wvfs))

    # === SG pipeline ===
    # savitzky golay filter: takes derivative of waveform plus smoothing
    sgflt_savitz = SavitzkyGolayFilter(sg_window_length, sg_flt_degree, 1)
    wvfs_sg = sgflt_savitz.(wvfs)

    # trigger finding on SG-filtered waveforms
    intflt_sg = IntersectMaximum(min_tot_intersect, max_tot_intersect)
    inters_thres = thresholdstats_mad.(wvfs_sg, min_threshold, max_threshold)
    inters = intflt_sg.(wvfs_sg, n_σ_threshold .* inters_thres)

    # integrate derivative
    integrator_filter = IntegratorFilter(gain=1)
    wvfs_int = integrator_filter.(wvfs_sg)

    # get blstats on integrated waveforms
    time_min = minimum(wvfs_int[1].time)
    time_max = last(wvfs_int[1].time)
    # baseline window ends at the first trigger (end of waveform if none), but spans at least 3 samples
    t_bl_min = time_min + 3*step(wvfs_int[1].time)
    t_bl_stop = map(x -> oftype(time_max, max(minimum(x; init = time_max), t_bl_min)), inters.x)
    bl_stats = signalstats.(wvfs_int, time_min, t_bl_stop)
    sigstats = signalstats.(wvfs_int, time_min, time_max)
    
    # flip around x-axis for discharge detection
    wvfs_int_flip = multiply_waveform.(wvfs_int, -1.0)
    inters_thres_DC = thresholdstats_mad.(wvfs_int_flip, min_dc_threshold, max_dc_threshold)
    inters_DC = intflt_sg.(wvfs_int_flip, n_σ_dc_threshold .* inters_thres_DC)

    # === Trap pipeline ===
    # PZ correction on integrated waveforms
    pz_filter = InvCRFilter(trap_pz_tau)
    wvfs_pz = pz_filter.(wvfs_int)

    # trapezoidal charge filter
    trap_filter = TrapezoidalChargeFilter(trap_rt, trap_ft)
    wvfs_trap = trap_filter.(wvfs_pz)

    # trigger finding on trap-filtered waveforms
    intflt_trap = IntersectMaximum(trap_min_tot_intersect, trap_max_tot_intersect)
    inters_thres_trap = thresholdstats_mad.(wvfs_trap, trap_min_threshold, trap_max_threshold)
    inters_trap = intflt_trap.(wvfs_trap, trap_n_σ_threshold .* inters_thres_trap)

    # discharge detection on flipped integrated waveforms 
    inters_thres_DC_trap = thresholdstats_mad.(wvfs_int_flip, trap_min_dc_threshold, trap_max_dc_threshold)
    inters_DC_trap = intflt_sg.(wvfs_int_flip, trap_n_σ_dc_threshold .* inters_thres_DC_trap)

    # === MWA pipeline ===
    # photon trace: the baseline-subtracted waveform through the inverse kernel of the template step; a channel without
    # template pars takes the derivative (kernel [1, -1]) instead
    t_axis = first(wvfs).time
    dt = step(t_axis)
    !has_kernel || pars.kernel.dt ≈ dt || throw(ArgumentError("The kernel taps are sampled at $(pars.kernel.dt), the waveforms at $dt"))
    photon_filter = CenteredConvolutionFilter(has_kernel ? Float64.(pars.kernel.taps) : [0.0, 1.0, -1.0])
    wvfs_bl = shift_waveform.(wvfs, [-median(wf.signal) for wf in wvfs])
    wvfs_photons = photon_filter.(wvfs_bl)

    # short (trigger) and long (energy) shaping, each scaled so that a unit stick (one photon after the kernel) has peak 1
    mwa_trigger_filter = MovingWindowMultiFilter(mwa_wl_trigger; passes = mwa_passes)
    mwa_energy_filter  = MovingWindowMultiFilter(mwa_wl_energy; passes = mwa_passes)
    stick = RDWaveform(t_axis, [i == length(t_axis) ÷ 2 ? 1.0 : 0.0 for i in eachindex(t_axis)])
    response_trigger = mwa_trigger_filter(stick).signal
    response_energy  = mwa_energy_filter(stick).signal
    wvfs_trigger = multiply_waveform.(mwa_trigger_filter.(wvfs_photons), 1 / maximum(response_trigger))
    wvfs_energy  = multiply_waveform.(mwa_energy_filter.(wvfs_photons), 1 / maximum(response_energy))
    wvfs_trigger = shift_waveform.(wvfs_trigger, [-median(wf.signal) for wf in wvfs_trigger])
    wvfs_energy  = shift_waveform.(wvfs_energy, [-median(wf.signal) for wf in wvfs_energy])

    # sites: the extrema beyond n σ of the trigger trace (matrix method) or of the energy trace
    wvfs_sites = mwa_solve_overlap ? wvfs_trigger : wvfs_energy
    mwa_thresholds = thresholdstats_mad.(wvfs_sites, mwa_min_threshold, mwa_max_threshold)
    sites = extrema_over_threshold.(wvfs_sites, mwa_n_σ_threshold .* mwa_thresholds)

    # matrix method: the amplitudes of all sites solved on the long shaping (otherwise the extrema of the energy trace)
    if mwa_solve_overlap
        solver = SiteSolver(response_energy ./ maximum(response_energy), argmax(response_energy), mwa_merge_distance * uconvert(NoUnits, mwa_wl_energy / dt))
        sites = map(zip(wvfs_energy, sites)) do (wf, s)
            sol = solver(wf.signal, uconvert.(NoUnits, (s.x .- first(t_axis)) ./ dt) .+ 1, uconvert.(NoUnits, s.x_tot ./ dt), Int.(sign.(s.max)))
            (x = convert.(eltype(t_axis), first(t_axis) .+ (sol.pos .- 1) .* dt), x_tot = convert.(eltype(t_axis), sol.tot .* dt), max = sol.A)
        end
    end

    # without the shaping transients at both ends of the trace
    mwa_margin = 2 * mwa_passes * max(mwa_wl_trigger, mwa_wl_energy)
    sites = map(sites) do s
        k = findall(t -> first(t_axis) + mwa_margin <= t <= last(t_axis) - mwa_margin, s.x)
        (x = s.x[k], x_tot = s.x_tot[k], max = s.max[k])
    end

    # output Table 
    return TypedTables.Table(
        blfc = blfc, timestamp = ts, eventID_fadc = evID, e_fc = efc,
        t_max = uconvert.(u"µs", estats.tmax), t_min = uconvert.(u"µs", estats.tmin), t_max_lar = uconvert.(u"µs", estats_trunc.tmax), t_min_lar = uconvert.(u"µs", estats_trunc.tmin),
        e_max = estats.max, e_min = estats.min, e_max_lar = estats_trunc.max, e_min_lar = estats_trunc.min,
        blmean = bl_stats.mean, blsigma = bl_stats.sigma, blslope = bl_stats.slope, bloffset = bl_stats.offset, 
        wfmean = sigstats.mean, wfsigma = sigstats.sigma, wfslope = sigstats.slope, wfoffset = sigstats.offset,
        # SG triggers
        threshold_sg = inters_thres, threshold_DC_sg = inters_thres_DC,
        trig_pos_sg = VectorOfVectors(inters.x), trig_max_sg = VectorOfVectors(inters.max),
        trig_pos_DC_sg = VectorOfVectors(inters_DC.x), trig_max_DC_sg = VectorOfVectors(inters_DC.max),
        # Trap triggers
        threshold_trap = inters_thres_trap, threshold_DC_trap = inters_thres_DC_trap,
        trig_pos_trap = VectorOfVectors(inters_trap.x), trig_pos_high_trap = VectorOfVectors(inters_trap.x_high),
        trig_pos_tot_trap = VectorOfVectors(inters_trap.x_tot), trig_max_trap = VectorOfVectors(inters_trap.max),
        trig_pos_DC_trap = VectorOfVectors(inters_DC_trap.x), trig_pos_high_DC_trap = VectorOfVectors(inters_DC_trap.x_high),
        trig_pos_tot_DC_trap = VectorOfVectors(inters_DC_trap.x_tot), trig_max_DC_trap = VectorOfVectors(inters_DC_trap.max),
        # MWA triggers
        threshold_mwa = mwa_thresholds, trig_pos_mwa = VectorOfVectors([s.x for s in sites]),
        trig_pos_tot_mwa = VectorOfVectors([s.x_tot for s in sites]), trig_max_mwa = VectorOfVectors([s.max for s in sites])
    )
end
export dsp_sipm_compressed
