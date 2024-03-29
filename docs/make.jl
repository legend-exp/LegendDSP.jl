# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using LegendDSP

# Doctest setup
DocMeta.setdocmeta!(
    LegendDSP,
    :DocTestSetup,
    :(using LegendDSP);
    recursive=true,
)

makedocs(
    sitename = "LegendDSP",
    modules = [LegendDSP],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://legend-exp.github.io/LegendDSP.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
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
