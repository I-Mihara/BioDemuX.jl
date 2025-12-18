using BioDemuX
using Documenter

DocMeta.setdocmeta!(BioDemuX, :DocTestSetup, :(using BioDemuX); recursive=true)

makedocs(;
    modules=[BioDemuX],
    authors="I.Mihara <issei.mihara@xforestx.com> and contributors",
    sitename="BioDemuX.jl",
    format=Documenter.HTML(;
        canonical="https://I-Mihara.github.io/BioDemuX.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Usage" => "usage.md",
        "CLI" => "cli.md",
        "Options" => "options.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/I-Mihara/BioDemuX.jl",
    devbranch="main",
)
