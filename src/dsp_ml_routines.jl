"""
    get_qc_classifier(wvfs::ArrayOfRDWaveforms, f_evaluate_qc::Function)
    get_qc_classifier(wvfs::ArrayOfRDWaveforms, f_evaluate_qc::Function, config::DSPConfig)

Quality-cut labels of the waveforms `wvfs`, evaluated by the support vector machine function `f_evaluate_qc`.

Each waveform is filtered five times with `HaarAveragingFilter(2)`, which reduces the number of samples by a factor of 32, and divided by the larger of the absolute values of its minimum and maximum (waveforms that are zero everywhere are not scaled). The filtered signals form the columns of a matrix that is passed to `f_evaluate_qc`. The method with `config` subtracts the baseline mean in `config.bl_window` first; the other method uses `wvfs` as given.

# Arguments
- `f_evaluate_qc`: function that returns the tuple `(labels, decision_values)` for a matrix with one column per waveform, as returned by [`get_qc_ml_func`](@ref)

# Returns
- vector of labels, one per waveform
"""
function get_qc_classifier end
export get_qc_classifier

function get_qc_classifier(wvfs::ArrayOfRDWaveforms, f_evaluate_qc::Function)
    # create Haar filter
    haar_flt = HaarAveragingFilter(2)

    # filter waveforms with Haar filter 5 times
    wvfs_flt_haar5 = wvfs .|> haar_flt .|> haar_flt .|> haar_flt .|> haar_flt .|> haar_flt

    # normalize with max of absolute extrema values
    norm_fact = map(x -> max(abs(first(x)), abs(last(x))), extrema.(wvfs_flt_haar5.signal))
    replace!(norm_fact, 0.0 => one(first(norm_fact)))

    wvfs_flt_haar5 = multiply_waveform.(wvfs_flt_haar5, 1 ./ norm_fact)

    y_pred, _ = f_evaluate_qc(flatview(convert(VectorOfSimilarArrays, wvfs_flt_haar5.signal)))
    y_pred
end

function get_qc_classifier(wvfs::ArrayOfRDWaveforms, f_evaluate_qc::Function, config::DSPConfig)
    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(config.bl_window), rightendpoint(config.bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    get_qc_classifier(wvfs, f_evaluate_qc)
end

"""
    get_qc_classifier_compressed(wvfs::ArrayOfRDWaveforms, f_evaluate_qc::Function)
    get_qc_classifier_compressed(wvfs::ArrayOfRDWaveforms, f_evaluate_qc::Function, config::DSPConfig)

Variant of [`get_qc_classifier`](@ref) for presummed waveforms: each waveform is filtered twice instead of five times with `HaarAveragingFilter(2)`, which reduces the number of samples by a factor of 4.
"""
function get_qc_classifier_compressed end 
export get_qc_classifier_compressed

function get_qc_classifier_compressed(wvfs::ArrayOfRDWaveforms, f_evaluate_qc::Function)
    # create Haar filter
    haar_flt = HaarAveragingFilter(2)

    # filter waveforms with Haar filter 2 times
    wvfs_flt_haar2 = wvfs .|> haar_flt .|> haar_flt

    # normalize with max of absolute extrema values
    norm_fact = map(x -> max(abs(first(x)), abs(last(x))), extrema.(wvfs_flt_haar2.signal))
    replace!(norm_fact, 0.0 => one(first(norm_fact)))

    wvfs_flt_haar2 = multiply_waveform.(wvfs_flt_haar2, 1 ./ norm_fact)

    y_pred, _ = f_evaluate_qc(flatview(convert(VectorOfSimilarArrays, wvfs_flt_haar2.signal)))
    y_pred
end

function get_qc_classifier_compressed(wvfs::ArrayOfRDWaveforms, f_evaluate_qc::Function, config::DSPConfig)
    # get baseline mean, std and slope
    bl_stats = signalstats.(wvfs, leftendpoint(config.bl_window), rightendpoint(config.bl_window))

    # substract baseline from waveforms
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)

    get_qc_classifier_compressed(wvfs, f_evaluate_qc)
end
