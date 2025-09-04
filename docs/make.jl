using EndoBeams
using Documenter

DocMeta.setdocmeta!(EndoBeams, :DocTestSetup, :(using EndoBeams); recursive=true)

makedocs(;
    modules=[EndoBeams],
    authors="Beatrice Bisighini",
    # repo="https://github.com/bisighinibeatrice/EndoBeams.jl.git",
    sitename="EndoBeams.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bisighinibeatrice.github.io/EndoBeams.jl",
        assets=String[],
        disable_git=true
    ),
    pages=[
        "Home" => "index.md",
    ],
)
