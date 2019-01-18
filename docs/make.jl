using Documenter, Pseudospectra
include("../src/PseudospectraPlots.jl")

makedocs(modules = [Pseudospectra, PseudospectraPlots],
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == true),
         sitename = "Pseudospectra.jl",
         pages = [
             "Home" => "index.md",
             "Usage" => "usage.md",
             "Library" => Any[
                 "Public" => "lib/public.md",
                 "Internals" => "lib/internals.md",
                 "Example matrix generators" => "lib/demos.md"
             ]
         ]
         )

deploydocs(repo = "github.com/RalphAS/Pseudospectra.jl.git")
