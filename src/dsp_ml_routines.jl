
"""
    get_qc_columns(wvfs::ArrayOfRDWaveforms, f_evaluate_qc) -> NamedTuple

Output columns of the ML QC evaluator `f_evaluate_qc` for `wvfs`.

`f_evaluate_qc(wvfs)` takes the baseline-subtracted waveforms, extracts its own
features and returns a `NamedTuple` of vectors with one element per waveform.
It must include the column `qc_label` (integer labels); further columns, e.g.
the probability of every QC category, are passed on into the DSP output.
With `f_evaluate_qc = missing` the result is `(qc_label = fill(-1, length(wvfs)),)`.
"""
function get_qc_columns end
export get_qc_columns

function get_qc_columns(wvfs::ArrayOfRDWaveforms, f_evaluate_qc)
    cols = f_evaluate_qc(wvfs)
    cols isa NamedTuple ||
        throw(ArgumentError("the ML QC evaluator must return a NamedTuple of columns, got $(typeof(cols))"))
    haskey(cols, :qc_label) || throw(ArgumentError("the ML QC evaluator returned no qc_label column: $(keys(cols))"))
    for (name, col) in pairs(cols)
        length(col) == length(wvfs) ||
            throw(DimensionMismatch("the ML QC evaluator column $name has $(length(col)) rows for $(length(wvfs)) waveforms"))
    end
    return merge(cols, (qc_label = Int.(cols.qc_label),))
end
get_qc_columns(wvfs::ArrayOfRDWaveforms, ::Missing) = (qc_label = fill(-1, length(wvfs)),)

"""
    get_qc_classifier(wvfs::ArrayOfRDWaveforms, f_evaluate_qc)
    get_qc_classifier(wvfs::ArrayOfRDWaveforms, f_evaluate_qc, config::DSPConfig)

QC labels of `wvfs`: the `qc_label` column of [`get_qc_columns`](@ref). With
`config`, the mean over `config.bl_window` is subtracted from the waveforms first.
"""
function get_qc_classifier end
export get_qc_classifier

get_qc_classifier(wvfs::ArrayOfRDWaveforms, f_evaluate_qc) = get_qc_columns(wvfs, f_evaluate_qc).qc_label

function get_qc_classifier(wvfs::ArrayOfRDWaveforms, f_evaluate_qc, config::DSPConfig)
    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(config.bl_window), rightendpoint(config.bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    get_qc_classifier(wvfs, f_evaluate_qc)
end
