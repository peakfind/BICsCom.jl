using BICsCom
using Documenter

DocMeta.setdocmeta!(BICsCom, :DocTestSetup, :(using BICsCom); recursive=true)

makedocs(;
    modules=[BICsCom],
    authors="Jiayi Zhang",
    sitename="BICsCom.jl",
    format=Documenter.HTML(;
        canonical="https://peakfind.github.io/BICsCom.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/peakfind/BICsCom.jl",
    devbranch="main",
)
