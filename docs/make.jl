using BICsCom
using Documenter

# for CMakieExt
using CairoMakie

DocMeta.setdocmeta!(BICsCom, :DocTestSetup, :(using BICsCom); recursive=true)

const cmakieExt = Base.get_extension(BICsCom, :CMakieExt)

makedocs(;
    modules=[BICsCom, cmakieExt],
    authors="Jiayi Zhang",
    sitename="BICsCom.jl",
    format=Documenter.HTML(;
        canonical="https://peakfind.github.io/BICsCom.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/peakfind/BICsCom.jl",
    devbranch="main",
)
