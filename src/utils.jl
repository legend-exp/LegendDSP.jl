# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

"""
    _create_dsp_config(dsp_metadata::PropDicts.PropDict)

Create a `DSPConfig` from a `PropDict` of DSP metadata.

# Arguments
- `dsp_metadata::PropDicts.PropDict`: DSP metadata

# Returns
- `dsp_config::DSPConfig`: DSP configuration
"""
function _create_dsp_config(dsp_metadata::PropDicts.PropDict)
    return DSPConfig{Float64}(
        # pick-off time for ENC noise calculations
        dsp_metadata.enc_pickoff_trap,
        dsp_metadata.enc_pickoff_zac,
        dsp_metadata.enc_pickoff_cusp,
        
        # filter lengths for CUSP and ZAC filters
        dsp_metadata.flt_length_cusp,
        dsp_metadata.flt_length_zac,
        
        # ADC threshold for t0 determination
        dsp_metadata.t0_threshold,
        # in-trace pile-up rejector threshold in sigmas
        dsp_metadata.inTraceCut_std_threshold,
        # Savitzky Golay polynominal order for current extraction
        dsp_metadata.sg_flt_degree,


        # fit window for basline extraction
        dsp_metadata.bl_window.min .. dsp_metadata.bl_window.max,
        # fit window for decay time extraction
        dsp_metadata.tail_window.min .. dsp_metadata.tail_window.max,
        # fit window for current extraction
        dsp_metadata.current_window.min .. dsp_metadata.current_window.max,
        # Integration length for QDrift extraction
        (first(dsp_metadata.qdrift_int_length):1.0*unit(first(dsp_metadata.qdrift_int_length)):last(dsp_metadata.qdrift_int_length)),
        # Integration length for LQ extraction
        (first(dsp_metadata.lq_int_length):1.0*unit(first(dsp_metadata.qdrift_int_length)):last(dsp_metadata.lq_int_length)),

        # rise and flat-top time grid scan ranges for trapezoidal filter
        dsp_metadata.e_grid_trap.rt.start:dsp_metadata.e_grid_trap.rt.step:dsp_metadata.e_grid_trap.rt.stop,
        dsp_metadata.e_grid_trap.ft.start:dsp_metadata.e_grid_trap.ft.step:dsp_metadata.e_grid_trap.ft.stop,
        # rise and flat-top time grid scan ranges for ZAC filter
        dsp_metadata.e_grid_zac.rt.start:dsp_metadata.e_grid_zac.rt.step:dsp_metadata.e_grid_zac.rt.stop,
        dsp_metadata.e_grid_zac.ft.start:dsp_metadata.e_grid_zac.ft.step:dsp_metadata.e_grid_zac.ft.stop,
        # rise and flat-top time grid scan ranges for CUSP filter
        dsp_metadata.e_grid_cusp.rt.start:dsp_metadata.e_grid_cusp.rt.step:dsp_metadata.e_grid_cusp.rt.stop,
        dsp_metadata.e_grid_cusp.ft.start:dsp_metadata.e_grid_cusp.ft.step:dsp_metadata.e_grid_cusp.ft.stop,
        # window length grid scan range for SG filter in current determination
        dsp_metadata.a_grid_wl_sg.start:dsp_metadata.a_grid_wl_sg.step:dsp_metadata.a_grid_wl_sg.stop,
        
        # flt defaults
        dsp_metadata.flt_defaults,
        
        # flt parameters
        dsp_metadata.kwargs_pars)
end

function get_fltpars(pd::PropDict, flt::Symbol, dspconfig::DSPConfig)
    if flt == :sg
        return get(get(pd, flt, PropDict()), :wl, dspconfig.default_flt_param[flt])
    else
        if !(haskey(pd, flt))
            return dspconfig.default_flt_param[flt].rt, dspconfig.default_flt_param[flt].ft
        else
            return get(pd[flt], :rt, dspconfig.default_flt_param[flt].rt), get(pd[flt], :ft, dspconfig.default_flt_param[flt].ft)
        end
    end
end