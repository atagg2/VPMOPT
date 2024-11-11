using VPMOPT
using Documenter

DocMeta.setdocmeta!(VPMOPT, :DocTestSetup, :(using VPMOPT); recursive=true)

makedocs(;
    modules=[VPMOPT],
    authors="Andrew Tagg <andrew.tagg57@gmail.com> and contributors",
    repo="https://github.com/byuflowlab/VPMOPT.jl/blob/{commit}{path}#{line}",
    sitename="VPMOPT.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://byuflowlab.github.io/VPMOPT.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/byuflowlab/VPMOPT.jl",
    devbranch="main",
)
