using MacaulayMatrix
using Documenter
import Literate

const _TUTORIAL_DIR = joinpath(@__DIR__, "src", "tutorials")
const _OUTPUT_DIR   = joinpath(@__DIR__, "src", "generated")

function literate_directory()
    for filename in readdir(_TUTORIAL_DIR)
        path = joinpath(_TUTORIAL_DIR, filename)
        Literate.markdown(path, _OUTPUT_DIR; documenter = true)
        Literate.notebook(path, _OUTPUT_DIR; documenter = true)
        Literate.script(path, _OUTPUT_DIR; documenter = true)
    end
end

literate_directory()

DocMeta.setdocmeta!(MacaulayMatrix, :DocTestSetup, :(using MacaulayMatrix); recursive=true)

makedocs(;
    modules=[MacaulayMatrix],
    authors="Benoît Legat <benoit.legat@gmail.com> and contributors",
    repo="https://gitlab.esat.kuleuven.be/benoit.legat/MacaulayMatrix.jl/blob/{commit}{path}#{line}",
    sitename="MacaulayMatrix.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://gitlab.esat.kuleuven.be/benoit.legat/MacaulayMatrix.jl",
        edit_link="main",
        assets=String[],
    ),
    pages = [
        "Index" => "index.md",
        "Tutorials" => map(
            file -> joinpath("generated", file),
            filter(
                file -> endswith(file, ".md"),
                sort(readdir(_OUTPUT_DIR)),
            ),
        ),
    ],
)

deploydocs(
    repo   = "gitlab.esat.kuleuven.be:benoit.legat/MacaulayMatrix.jl.git",
    push_preview = true,
)
