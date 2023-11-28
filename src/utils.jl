"""
    create_dsp_config(dsp_metadata::PropDicts.PropDict)

Create a `DSPConfig` from a `PropDict` of DSP metadata.

# Arguments
- `dsp_metadata::PropDicts.PropDict`: DSP metadata

# Returns
- `dsp_config::DSPConfig`: DSP configuration
"""
function create_dsp_config end
export create_dsp_config


function create_dsp_config(dsp_metadata::PropDicts.PropDict)
    return DSPConfig{Float64}(
        # pick-off time for ENC noise calculations
        dsp_metadata.enc_pickoff_trap*u"µs",
        dsp_metadata.enc_pickoff_zac*u"µs",
        dsp_metadata.enc_pickoff_cusp*u"µs",
        # filter lengths for CUSP and ZAC filters
        dsp_metadata.flt_length_cusp*u"µs",
        dsp_metadata.flt_length_zac*u"µs",
        # in-trace pile-up rejector threshold in sigmas
        dsp_metadata.inTraceCut_std_threshold,
        # fit window for basline extraction
        (dsp_metadata.bl_mean.min*u"µs", dsp_metadata.bl_mean.max*u"µs"),
        # fit window for decay time extraction
        (dsp_metadata.pz_fit.min*u"µs", dsp_metadata.pz_fit.max*u"µs"),
        # ADC threshold for t0 determination
        dsp_metadata.t0_threshold,
        # rise and flat-top time grid scan ranges for trapezoidal filter
        dsp_metadata.e_grid_trap.rt.start*u"µs":dsp_metadata.e_grid_trap.rt.step*u"µs":dsp_metadata.e_grid_trap.rt.stop*u"µs",
        dsp_metadata.e_grid_trap.ft.start*u"µs":dsp_metadata.e_grid_trap.ft.step*u"µs":dsp_metadata.e_grid_trap.ft.stop*u"µs",
        # rise and flat-top time grid scan ranges for ZAC filter
        dsp_metadata.e_grid_zac.rt.start*u"µs":dsp_metadata.e_grid_zac.rt.step*u"µs":dsp_metadata.e_grid_zac.rt.stop*u"µs",
        dsp_metadata.e_grid_zac.ft.start*u"µs":dsp_metadata.e_grid_zac.ft.step*u"µs":dsp_metadata.e_grid_zac.ft.stop*u"µs",
        # rise and flat-top time grid scan ranges for CUSP filter
        dsp_metadata.e_grid_cusp.rt.start*u"µs":dsp_metadata.e_grid_cusp.rt.step*u"µs":dsp_metadata.e_grid_cusp.rt.stop*u"µs",
        dsp_metadata.e_grid_cusp.ft.start*u"µs":dsp_metadata.e_grid_cusp.ft.step*u"µs":dsp_metadata.e_grid_cusp.ft.stop*u"µs",
        # window length grid scan range for SG filter in current determination
        dsp_metadata.a_grid_wl_sg.start*u"ns":dsp_metadata.a_grid_wl_sg.step*u"ns":dsp_metadata.a_grid_wl_sg.stop*u"ns")
end
