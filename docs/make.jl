using Documenter, BioSequenceGraphs

makedocs(
    format = Documenter.HTML(),
    sitename = "GenomeGraphs.jl",
    authors = "Ben J. Ward & Arda Akdemir",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Guide" => "man/guide.md"
        ]
    ],
    
)