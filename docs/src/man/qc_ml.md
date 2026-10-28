# Machine-Learning Quality Cuts

Some recorded waveforms do not show a clean pulse, for example because of electrical discharges.
The DSP chains can attach to every waveform a label from a classifier that was trained on waveforms of known class.
The label is stored in the column `qc_label` of the chain output, and cuts on this column and on other DSP parameters are made later in the analysis.
This page describes how the classifier is trained, how it is evaluated on waveforms, and how it enters the DSP chains.

## Overview

The classifier is a support vector machine (SVM) of the LIBSVM library with a radial basis function kernel.
LegendDSP.jl provides three functions for it:

- [`get_qc_ml_func`](@ref) trains the SVM on a matrix of features and a vector of labels and returns the evaluation function, called `f_evaluate_qc` in the rest of the package.
- [`get_qc_classifier`](@ref) converts waveforms to features and calls `f_evaluate_qc`. It is meant for waveforms at the full sampling rate.
- [`get_qc_classifier_compressed`](@ref) does the same for the presummed waveforms of the compressed data.

## Features

[`get_qc_classifier`](@ref) and [`get_qc_classifier_compressed`](@ref) use the following features, which have to be the same as those of the training set:

1. The waveform is reduced in length by applying [`HaarAveragingFilter`](@ref)`(2)` repeatedly: five times for `get_qc_classifier` and two times for `get_qc_classifier_compressed`. Each application halves the number of samples. Waveforms of 4096 samples have 128 samples after five applications.
2. The waveform is divided by the larger of the absolute values of its minimum and maximum. A waveform for which this value is zero is divided by 1.
3. The normalized samples of one waveform form one column of the feature matrix, so that the matrix has one row per sample and one column per waveform.

The baseline of the waveforms has to be subtracted first, and the methods that take a [`DSPConfig`](@ref) as a third argument do this: `get_qc_classifier(wvfs, f_evaluate_qc, config)` subtracts the mean of each waveform in `config.bl_window` and then calls `get_qc_classifier(wvfs, f_evaluate_qc)`.
All waveforms of one call need the same number of samples.

## The Evaluation Function

`f_evaluate_qc` is any function that takes the feature matrix and returns a tuple `(labels, decision_values)`; only the first element is used.
The function that [`get_qc_ml_func`](@ref) returns is `Base.Fix1(svmpredict, model)` of LIBSVM, which has this form.
[`get_qc_classifier`](@ref) and [`get_qc_classifier_compressed`](@ref) return the vector of predicted labels.
The DSP chains convert the labels to integers with `Int.(...)`, so the labels must be integer-valued.

## Training a Classifier

[`get_qc_ml_func`](@ref)`(dwts_norm, dc_labels, hyperparams)` takes the training features `dwts_norm` as a matrix with one column per training waveform, the labels `dc_labels` as a vector with one label per column, and the hyperparameters as a `PropDict` with the following keys:

| Key | Meaning |
|:---|:---|
| `cost` | cost parameter of the SVM |
| `gamma` | width parameter of the radial basis function kernel |
| `weights` | `PropDict` that maps each label, written as a key such as `Symbol("0")`, to the weight of that class |
| `probability` | whether to train a model for probability estimates |
| `cache_size` | size of the kernel cache in MB |
| `coef0` | kernel constant; not used by the radial basis function kernel but passed on to LIBSVM |
| `shrinking` | whether to use the shrinking heuristic |
| `tolerance` | tolerance of the termination criterion |

The following example trains a classifier on synthetic waveforms of two classes: pulses without and with a superimposed damped oscillation after the pulse.
The features of the training set are obtained with the same function that is used for the evaluation.
It is called with a function that stores the feature matrix instead of classifying.

```@example qc_ml
using LegendDSP
using PropDicts
using RadiationDetectorDSP
using RadiationDetectorSignals
using StableRNGs
using Unitful

function pulse(rng, t; t0, amplitude, glitch, rise = 0.3, tau = 50.0, offset = 10.0, noise = 1.0)
    x = ustrip.(u"µs", t)
    signal = [xi < t0 ? 0.0 : amplitude * (1 - exp(-(xi - t0) / rise)) * exp(-(xi - t0) / tau) for xi in x]
    ripple = [xi < t0 + 5 ? 0.0 : sin(2π * xi / 3.0) * exp(-(xi - t0 - 5) / 20) for xi in x]
    RDWaveform(t, offset .+ signal .+ glitch .* ripple .+ noise .* randn(rng, length(x)))
end

rng = StableRNG(11)
t = range(0.0u"µs", step = 0.016u"µs", length = 4096)

function waveforms(n; glitch)
    ArrayOfRDWaveforms([
        pulse(rng, t; t0 = 20 + 10 * rand(rng), amplitude = 50 + 100 * rand(rng), glitch = glitch * (0.5 + rand(rng)))
        for _ in 1:n
    ])
end

function subtract_baseline(wvfs)
    bl_stats = signalstats.(wvfs, 0u"µs", 15u"µs")
    shift_waveform.(wvfs, -bl_stats.mean)
end

wvfs_train = ArrayOfRDWaveforms(collect(vcat(waveforms(100; glitch = 0.0), waveforms(100; glitch = 30.0))))
labels_train = vcat(fill(0.0, 100), fill(1.0, 100))

features = Ref(zeros(0, 0))
record_features(X) = (features[] = Matrix{Float64}(X); (zeros(size(X, 2)), zeros(size(X, 2))))
get_qc_classifier(subtract_baseline(wvfs_train), record_features)
size(features[])
```

The feature matrix has 128 rows, one for each sample after five Haar averaging steps, and 200 columns, one for each training waveform.

```@example qc_ml
hyperparams = PropDict(
    :cost => 1.0, :gamma => 0.1, :weights => PropDict(Symbol("0") => 1.0, Symbol("1") => 1.0),
    :probability => false, :cache_size => 200, :coef0 => 0.0, :shrinking => true, :tolerance => 0.001,
)
f_evaluate_qc = get_qc_ml_func(features[], labels_train, hyperparams)
```

The classifier is applied to a new set of waveforms, 20 of each class:

```@example qc_ml
wvfs_test = ArrayOfRDWaveforms(collect(vcat(waveforms(20; glitch = 0.0), waveforms(20; glitch = 30.0))))
labels_test = get_qc_classifier(subtract_baseline(wvfs_test), f_evaluate_qc)
(count(==(0.0), labels_test[1:20]), count(==(1.0), labels_test[21:40]))
```

## Use in the DSP Chains

The functions that compute `qc_label` take the evaluation function as the keyword argument `f_evaluate_qc`, or as the last positional argument for [`dsp_qc_flt_optimization`](@ref) and [`dsp_qc_flt_optimization_compressed`](@ref).
Its type is `Union{Function, Missing}`:

- [`dsp_icpc`](@ref), [`dsp_sg_optimization`](@ref), and [`dsp_mwa_optimization`](@ref) call [`get_qc_classifier`](@ref).
- [`dsp_icpc_compressed`](@ref), [`dsp_sg_optimization_compressed`](@ref), and [`dsp_mwa_optimization_compressed`](@ref) call [`get_qc_classifier_compressed`](@ref) on the presummed waveforms.
- [`dsp_qc_flt_optimization`](@ref) and [`dsp_qc_flt_optimization_compressed`](@ref) call one of the two functions as well, on waveforms to which the baseline subtraction and the pole-zero correction have been applied. The other functions call it on baseline-subtracted waveforms without pole-zero correction.
- If `f_evaluate_qc` is `missing`, no classification is made, and `qc_label` is `-1` for all events.

```julia
f_evaluate_qc = get_qc_ml_func(features, labels, hyperparams)

dsp_out = dsp_icpc_compressed(data, config, τ, pars_filter; f_evaluate_qc)
dsp_out.qc_label
```

Juleana.jl creates the evaluation function in a single place, where it is set to `missing`, so that the ML step is skipped and `qc_label` is `-1` in the production output (see [Use in Juleana.jl](@ref)).

The functions `dsp_icpc_compressed` and the other chains are described in [HPGe DSP Chains](@ref), and the use of `qc_label` in the filter scans in [Filter Optimization](@ref).
