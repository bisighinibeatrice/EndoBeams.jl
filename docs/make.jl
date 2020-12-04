using EndoBeams
using Documenter

makedocs(;
    modules=[EndoBeams],
    authors="Beatrice Bisighini, Baptiste Pierrat, Miquel Aguirre",
    repo="https://gitlab.emse.fr/pierrat/EndoBeams.jl/blob/{commit}{path}#L{line}",
    sitename="EndoBeams.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="http://pierrat.gitlab.emse.fr/EndoBeams.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
