using Documenter, Pseudospectra
include("../src/PseudospectraPlots.jl")

makedocs(modules = [Pseudospectra, PseudospectraPlots],
         format = :html,
         sitename = "Pseudospectra",
         pages = Any["Home" => "index.md",
                     # "Usage" => ...
                     "Library" => "library.md",
                     "Internals" => "internals.md"
                     ]
         )
