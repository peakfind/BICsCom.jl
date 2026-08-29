using BICsCom
using Documenter, DocumenterCitations

# for CMakieExt
using CairoMakie

DocMeta.setdocmeta!(BICsCom, :DocTestSetup, :(using BICsCom); recursive=true)

const CMakieExt = Base.get_extension(BICsCom, :CMakieExt)

bib = CitationBibliography("src/refs.bib")

makedocs(;
    modules=[BICsCom, CMakieExt],
    authors="Jiayi Zhang",
    sitename="BICsCom.jl",
    format=Documenter.HTML(;
        canonical="https://peakfind.github.io/BICsCom.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
        "References" => "references.md"
    ],
    plugins=[bib],
)

deploydocs(;
    repo="github.com/peakfind/BICsCom.jl",
    devbranch="main",
)
