using Documenter
using EndoBeams   # replace with your module name

makedocs(
    sitename = "EndoBeams Documentation",
    modules = [EndoBeams],
    pages = [
        "Home" => "index.md"
    ]
)

deploydocs(
    repo = "https://github.com/bisighinibeatrice/EndoBeams.jl.git",
    branch = "gh-pages",
)