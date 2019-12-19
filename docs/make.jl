using Documenter, GenomeGraphs

makedocs(
    modules = [GenomeGraphs],
    format = Documenter.HTML(),
    sitename = "GenomeGraphs.jl",
    authors = "Ben J. Ward & Arda Akdemir",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Guide" => "man/guide.md"
        ],
        "API" => [
            "Graphs" => "api/SequenceDistanceGraph.md"
        ]
    ],
    
)

deploydocs(
    repo = "github.com/BioJulia/GenomeGraphs.jl.git",
    deps = nothing,
    make = nothing
)