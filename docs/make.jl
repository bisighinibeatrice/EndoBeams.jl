using EndoBeams
using Documenter

DocMeta.setdocmeta!(EndoBeams, :DocTestSetup, :(using EndoBeams); recursive=true)

makedocs(;
    modules=[EndoBeams],
    authors="Baptiste Pierrat",
    repo="https://gitlab.emse.fr/pierrat/EndoBeams.jl/blob/{commit}{path}#{line}",
    sitename="EndoBeams.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pierrat.gitlab.emse.fr/EndoBeams.jl",
        assets=String[],
        disable_git=true
    ),
    pages=[
        "Home" => "index.md",
    ],
)
