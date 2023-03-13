using Macaulay
using Documenter

DocMeta.setdocmeta!(Macaulay, :DocTestSetup, :(using Macaulay); recursive=true)

makedocs(;
    modules=[Macaulay],
    authors="Benoît Legat <benoit.legat@gmail.com> and contributors",
    repo="https://gitlab.esat.kuleuven.be/benoit.legat/Macaulay.jl/blob/{commit}{path}#{line}",
    sitename="Macaulay.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://gitlab.esat.kuleuven.be/benoit.legat/Macaulay.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
