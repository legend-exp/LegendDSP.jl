"""
    get_qc_ml_func(dwts_norm::Matrix{<:Real}, dc_labels::Vector{<:Real}, hyperparams::PropDict)

Create a function that takes a signal and returns the prediction using the given hyperparameters and discrete wavelet transforms.
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