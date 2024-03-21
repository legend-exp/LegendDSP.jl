"""
    get_qc_classifier(wvfs::ArrayOfRDWaveforms, f_evaluate_qc::Function)

Get a classifier for quality cuts using the given waveforms and evaluation function which evulates a predefined SVM model.
"""
function get_qc_classifier(wvfs::ArrayOfRDWaveforms, f_evaluate_qc::Function)
    # create Haar filter
    haar_flt = HaarAveragingFilter(2)

    # filter waveforms with Haar filter 5 times
    wvfs_flt_haar5 = wvfs .|> haar_flt .|> haar_flt .|> haar_flt .|> haar_flt .|> haar_flt

    # normalize with max of absolute extrema values
    norm_fact = map(x -> max(abs(first(x)), abs(last(x))), extrema.(wvfs_flt_haar5.signal))
    replace!(norm_fact, 0.0 => one(first(norm_fact)))

    wvfs_flt_haar5 = multiply_waveform.(wvfs_flt_haar5, 1 ./ norm_fact)

    y_pred, _ = f_evaluate_qc(flatview(VectorOfSimilarArrays(wvfs_flt_haar5.signal)))
    y_pred
end