# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using Literate
using LegendDSP

function fix_literate_output(content)
    content = replace(content, "EditURL = \"@__REPO_ROOT_URL__/\"" => "")
    return content
end

gen_content_dir = joinpath(@__DIR__, "src", "tutorials")
for tut_lit_fn in filter(fn -> endswith(fn, "_lit.jl"), readdir(gen_content_dir))
    lit_src_fn = joinpath(gen_content_dir, tut_lit_fn)
    tut_basename = tut_lit_fn[1:end-7] # remove "_lit.jl"
    Literate.markdown(lit_src_fn, gen_content_dir, name = tut_basename, documenter = true, credit = true, postprocess = fix_literate_output)
    Literate.notebook(lit_src_fn, gen_content_dir, execute = false, name = tut_basename, documenter = true, credit = true)
    Literate.script(lit_src_fn, gen_content_dir, keep_comments = false, name = tut_basename, documenter = true, credit = false)
end

# Doctest setup
DocMeta.setdocmeta!(
    LegendDSP,
    :DocTestSetup,
    :(using LegendDSP);
    recursive=true,
)

makedocs(
    sitename = "LegendDSP.jl",
    modules = [LegendDSP],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://legend-exp.github.io/LegendDSP.jl/stable/",
        size_threshold = nothing,
        size_threshold_warn = nothing,
        example_size_threshold = nothing,
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "Waveforms and Configuration" => "man/waveforms_config.md",
            "HPGe DSP Chains" => "man/dsp_chains.md",
            "Filters" => "man/filters.md",
            "Waveform Parameters" => "man/routines.md",
            "Decay Time" => "man/decay_time.md",
            "Filter Optimization" => "man/filter_optimization.md",
            "SiPM and PMT Processing" => "man/sipm_pmt.md",
            "Machine-Learning Quality Cuts" => "man/qc_ml.md",
            "Use in Juleana.jl" => "man/juleana.md",
        ],
        "Tutorials" => [
            "tutorials/energy.md",
            "tutorials/aoe.md",
        ],
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    warnonly = ("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/legend-exp/LegendDSP.jl.git",
    forcepush = true,
    push_preview = true,
)
