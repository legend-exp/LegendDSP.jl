"""
    get_qc_ml_func(dwts_norm::Matrix{<:Real}, dc_labels::Vector{<:Real}, hyperparams::PropDict)

Train a support vector machine with a radial basis function kernel and return a function that classifies waveforms.

# Arguments
- `dwts_norm`: training features, one column per waveform
- `dc_labels`: training label of each waveform
- `hyperparams`: `PropDict` with the `LIBSVM` settings `cost`, `gamma`, `weights` (class weights, keyed by the class label), `probability`, `cache_size`, `coef0`, `shrinking` and `tolerance`

# Returns
- function `f(x)` for a feature matrix `x` with one column per waveform. It returns the tuple `(labels, decision_values)` of `LIBSVM.svmpredict`. Pass it as `f_evaluate_qc` to the DSP chains and to [`get_qc_classifier`](@ref).
"""
function get_qc_ml_func(dwts_norm::Matrix{<:Real}, dc_labels::Vector{<:Real}, hyperparams::PropDict)
    # create SVM model from hyperparams and discrete wavelet transforms and dc_labels
    model = svmtrain(dwts_norm, dc_labels, 
                cost=hyperparams.cost, 
                kernel=LIBSVM.Kernel.RadialBasis, 
                gamma=hyperparams.gamma,
                weights = Dict(parse.(Float64, string.(keys(hyperparams.weights))) .=> values(hyperparams.weights)),
                probability=hyperparams.probability,
                cachesize=Float64(hyperparams.cache_size),
                coef0=hyperparams.coef0,
                shrinking=hyperparams.shrinking,
                tolerance=hyperparams.tolerance
                )

    # return a function that takes a signal and returns the prediction
    Base.Fix1(svmpredict, model)
end
export get_qc_ml_func