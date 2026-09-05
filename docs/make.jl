using BICsCom
using Documenter, DocumenterCitations, DocumenterCodeBlocks
using Literate

DocMeta.setdocmeta!(BICsCom, :DocTestSetup, :(using BICsCom); recursive=true)

bib = CitationBibliography("src/refs.bib")

function generate_tutorials(source_dir, output_dir; exclude=[])
    for file in readdir(source_dir)
        if endswith(file, ".jl") && !(file in exclude)
            input = joinpath(source_dir, file)
            output = output_dir
            Literate.markdown(input, output; documenter=true)
        end
    end
end

source = joinpath(@__DIR__, "src", "tutorials_literate")
output = joinpath(@__DIR__, "src", "tutorials")
generate_tutorials(source, output)

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
        "Introduction" => ["introduction/bics.md", 
                           "introduction/problem.md",
                           "introduction/method.md"],
        "Tutorials" => ["tutorials/standing_cylinder_te.md",
                        "tutorials/propagating_cylinder_te.md"],
        "API" => "api.md",
        "References" => "references.md"
    ],
    plugins=[bib, CodeBlocks()],
)

deploydocs(;
    repo="github.com/peakfind/BICsCom.jl",
    devbranch="main",
)
