using LegendDSP, PropDicts, Unitful, IntervalSets

function make_fake_config()
    pd = PropDict(

        :enc_pickoff_trap => 40.0u"µs",
        :enc_pickoff_zac  => 41.0u"µs",
        :enc_pickoff_cusp => 41.0u"µs",

        :bl_window => PropDict(
            :min => 0.0u"µs",
            :max => 39.0u"µs",
        ),
        :tail_window => PropDict(
            :min => 70.0u"µs",
            :max => 110.0u"µs",
        ),
        :current_window => PropDict(
            :min => 43.0u"µs",
            :max => 62.0u"µs",
        ),

        :auxpz1_window => PropDict(
            :min => 70.0u"µs",
            :max => 90.0u"µs",
        ),
        :auxpz2_window => PropDict(
            :min => 90.0u"µs",
            :max => 110.0u"µs",
        ),
        :flt_length_cusp => 38.0u"µs",
        :flt_length_zac  => 38.0u"µs",
        :t0_threshold             => 4.0,
        :inTraceCut_std_threshold => 5,
        :sg_flt_degree            => 3,
        :qdrift_int_length => (2.5u"µs":0.1u"µs":5.0u"µs"),
        :lq_int_length     => (2.5u"µs":0.1u"µs":5.0u"µs"),

        :e_grid_trap => PropDict(
            :rt => PropDict(
                :start => 1.0u"µs",
                :stop  => 16.0u"µs",
                :step  => 0.5u"µs",
            ),
            :ft => PropDict(
                :start => 1.0u"µs",
                :stop  => 4.0u"µs",
                :step  => 0.2u"µs",
            ),
        ),

        :e_grid_zac => PropDict(
            :rt => PropDict(
                :start => 1.0u"µs",
                :stop  => 16.0u"µs",
                :step  => 0.5u"µs",
            ),
            :ft => PropDict(
                :start => 1.0u"µs",
                :stop  => 4.0u"µs",
                :step  => 0.2u"µs",
            ),
        ),

        :e_grid_cusp => PropDict(
            :rt => PropDict(
                :start => 1.0u"µs",
                :stop  => 16.0u"µs",
                :step  => 0.5u"µs",
            ),
            :ft => PropDict(
                :start => 1.0u"µs",
                :stop  => 4.0u"µs",
                :step  => 0.2u"µs",
            ),
        ),

        :a_grid_wl_sg => PropDict(
            :start => 30.0u"ns",
            :stop  => 350.0u"ns",
            :step  => 32.0u"ns",
        ),

        :a_grid_wl_mwa => PropDict(
            :start => 16.0u"ns",
            :stop  => 576.0u"ns",
            :step  => 16.0u"ns",
        ),

        :flt_defaults => PropDict(
            :sg   => 100.0u"ns",
            :mwa  => 96.0u"ns",
            :trap => PropDict(:rt => 5.0u"µs", :ft => 2.5u"µs"),
            :zac  => PropDict(:rt => 5.0u"µs", :ft => 2.5u"µs"),
            :cusp => PropDict(:rt => 5.0u"µs", :ft => 2.5u"µs"),
        ),

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
    )

    return DSPConfig(pd)
end
