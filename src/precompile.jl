# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

# Not covered: `dsp_sipm_compressed` and `dsp_pmts` (they need bit-drop encoded waveforms and PMT metadata).
@setup_workload begin
    n_wvf = 2

    hpge_signal(n) = [i < 3000 ? 1000.0 : i < 3125 ? 1000.0 + 10000.0 * (i - 3000) / 125 :
                      1000.0 + 10000.0 * exp(-(i - 3125) / 31250) for i in 1:n]
    hpge_wvfs() = ArrayOfRDWaveforms([RDWaveform((0:8191) .* 16.0u"ns", hpge_signal(8192)) for _ in 1:n_wvf])

    sipm_signal(n) = [i < 3125 ? 0.0 : 5.0 * exp(-(i - 3125) / 30.0) for i in 1:n]
    sipm_wvfs() = ArrayOfRDWaveforms([RDWaveform((0:6249) .* 16.0u"ns", sipm_signal(6250)) for _ in 1:n_wvf])

    window(lo, hi) = PropDict(:min => lo, :max => hi)
    grid(start, stop, step) = PropDict(:start => start, :stop => stop, :step => step)
    eflt_grid = PropDict(:rt => grid(1.0u"µs", 16.0u"µs", 1.0u"µs"), :ft => grid(1.0u"µs", 4.0u"µs", 0.5u"µs"))
    flt_default = PropDict(:rt => 5.0u"µs", :ft => 2.5u"µs")

    config = DSPConfig(PropDict(
        :enc_pickoff_trap => 40.0u"µs",
        :enc_pickoff_zac  => 41.0u"µs",
        :enc_pickoff_cusp => 41.0u"µs",
        :bl_window        => window(0.0u"µs", 39.0u"µs"),
        :tail_window      => window(70.0u"µs", 110.0u"µs"),
        :current_window   => window(43.0u"µs", 62.0u"µs"),
        :auxpz1_window    => window(70.0u"µs", 90.0u"µs"),
        :auxpz2_window    => window(90.0u"µs", 110.0u"µs"),
        :flt_length_cusp  => 38.0u"µs",
        :flt_length_zac   => 38.0u"µs",
        :t0_threshold             => 4.0,
        :inTraceCut_std_threshold => 5,
        :sg_flt_degree            => 3,
        :qdrift_int_length => (2.5u"µs":0.5u"µs":5.0u"µs"),
        :lq_windows        => [[2.5, 0.0, 2.5]u"µs", [0.5, 1.0, 1.5]u"µs"],
        :e_grid_trap => eflt_grid,
        :e_grid_zac  => eflt_grid,
        :e_grid_cusp => eflt_grid,
        :a_grid_wl_sg  => grid(30.0u"ns", 350.0u"ns", 64.0u"ns"),
        :a_grid_wl_mwa => grid(16.0u"ns", 576.0u"ns", 112.0u"ns"),
        :a_mwa_windows => [48.0u"ns", 288.0u"ns"],
        :a_mwa_passes  => 3,
        :a_sg_windows  => [60.0u"ns", 100.0u"ns"],
        :a_wsg_windows => [112.0u"ns"],
        :a_wsg_weights => [2],
        :a_wsg_degree  => 2,
        :flt_defaults  => PropDict(:sg => 100.0u"ns", :mwa => 96.0u"ns",
                                   :trap => flt_default, :zac => flt_default, :cusp => flt_default, :dplms => flt_default),
        :kwargs_pars => PropDict(
            :fc_bit_depth             => 16,
            :t0_flt_pars              => [40.0u"ns", 100.0u"ns", 2000.0u"ns"],
            :t0_mintot                => 1500.0u"ns",
            :tx_mintot                => 32.0u"ns",
            :intrace_mintot           => 100.0u"ns",
            :int_interpolation_order  => 3,
            :int_interpolation_length => 100.0u"ns",
            :sig_interpolation_order  => 3,
            :sig_interpolation_length => 700.0u"ns",
        ),
    ))

    sipm_filter(extra...) = PropDict(
        :n_σ_threshold => 3.0, :min_threshold => -1.0, :max_threshold => 1.0,
        :n_σ_dc_threshold => 5.0, :min_dc_threshold => -4.0, :max_dc_threshold => 4.0,
        :min_tot_intersect => 70.0u"ns", :max_tot_intersect => 150.0u"ns", extra...)
    sipm_config = PropDict(
        :t0_hpge_window => [47.0u"µs", 53.0u"µs"],
        :sg_flt_degree  => 3,
        :filters => PropDict(:sg => sipm_filter(), :trap => sipm_filter(:rt => 100.0u"ns", :ft => 50.0u"ns", :pz_tau => 3.0u"µs"),
                          :xt_tag => sipm_filter(:window_length => 500.0u"ns")),
    )
    sipm_sg_config = sipm_filter(:sg_flt_degree => 3)
    sipm_pars = PropDict(:sg => PropDict(:wl => 200.0u"ns"))
    sipm_opt_config = PropDict(
        :e_grid_wl => [100.0u"ns", 200.0u"ns"],
        :threshold => PropDict(:min_cut => -1.0, :max_cut => 1.0, :n_wvfs => n_wvf),
    )

    τ = 500.0u"µs"
    pars_trap = PropDict(:trap => flt_default)
    blmean = fill(1000.0, n_wvf)

    @compile_workload begin
        windowed = hpge_wvfs()
        presummed = hpge_wvfs()
        baseline = fill(0.0f0, n_wvf)
        ts = fill(UInt64(0), n_wvf)
        evnum = UInt32.(1:n_wvf)
        daqe = fill(UInt16(0), n_wvf)

        hpge_data = Table(waveform = windowed, baseline = baseline, timestamp = ts, eventnumber = evnum, daqenergy = daqe)
        hpge_data_compressed = Table(
            waveform_windowed = windowed, waveform_presummed = presummed, presum_rate = fill(UInt16(1), n_wvf),
            baseline = baseline, timestamp = ts, eventnumber = evnum, daqenergy = daqe,
            t_sat_lo = fill(UInt16(0), n_wvf), t_sat_hi = fill(UInt16(0), n_wvf), deadtime = fill(UInt16(0), n_wvf),
        )
        pulser_data = Table(waveform = presummed, baseline = baseline, timestamp = ts, eventnumber = evnum, daqenergy = daqe)

        dsp_icpc(hpge_data, config, τ, PropDict())
        dsp_icpc_compressed(hpge_data_compressed, config, τ, PropDict())
        dsp_aux(hpge_data, config)
        dsp_aux_compressed(hpge_data_compressed, config)
        dsp_puls(pulser_data, config)
        dsp_puls_compressed(hpge_data_compressed, config)

        dsp_decay_times(windowed, config)
        dsp_qc_flt_optimization(windowed, config, τ, missing)
        dsp_qc_flt_optimization_compressed(windowed, config, τ, missing)
        dsp_qdrift_flt_optimization(windowed, blmean, config, τ)
        dsp_trap_rt_optimization(windowed, config, τ)
        dsp_cusp_rt_optimization(windowed, config, τ)
        dsp_zac_rt_optimization(windowed, config, τ)
        dsp_trap_ft_optimization(windowed, config, τ, 5.0u"µs")
        dsp_cusp_ft_optimization(windowed, config, τ, 5.0u"µs")
        dsp_zac_ft_optimization(windowed, config, τ, 5.0u"µs")
        dsp_sg_optimization(windowed, config, τ, pars_trap)
        dsp_sg_optimization_compressed(windowed, presummed, config, τ, pars_trap; presum_rate = 1.0)
        dsp_mwa_optimization(windowed, config, τ, pars_trap)
        dsp_mwa_optimization_compressed(windowed, presummed, config, τ, pars_trap; presum_rate = 1.0)

        sipm = sipm_wvfs()
        sipm_data = Table(waveform = sipm, baseline = baseline, timestamp = ts, eventnumber = evnum, daqenergy = daqe)
        dsp_sipm(sipm_data, sipm_config, sipm_pars)
        dsp_sg_sipm_thresholds_compressed(sipm, 200.0u"ns", sipm_config)
        dsp_sg_sipm_optimization_compressed(sipm, sipm_sg_config, sipm_opt_config)
    end
end
