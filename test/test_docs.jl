# This file is a part of LegendDSP.jl, licensed under the MIT License (MIT).

using Test
using LegendDSP
import Documenter

Documenter.DocMeta.setdocmeta!(
    LegendDSP,
    :DocTestSetup,
    :(using LegendDSP);
    recursive=true,
)
Documenter.doctest(LegendDSP)
