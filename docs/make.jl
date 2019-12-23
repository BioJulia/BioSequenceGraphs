using Documenter, GenomeGraphs

makedocs(
    modules = [GenomeGraphs, GenomeGraphs.Graphs, GenomeGraphs.MerTools],
    format = Documenter.HTML(),
    sitename = "GenomeGraphs.jl",
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Guide" => "man/guide.md"
        ],
        "API" => [
            "Graphs submodule" => "api/Graphs.md"
        ]
    ],
    
)

deploydocs(
    repo = "github.com/BioJulia/GenomeGraphs.jl.git",
    push_preview = true,
    deps = nothing,
    make = nothing
)