using Documenter, Pseudospectra
include("../src/PseudospectraPlots.jl")

makedocs(modules = [Pseudospectra, PseudospectraPlots],
         format = :html,
         sitename = "Pseudospectra.jl",
         pages = Any[
             "Home" => "index.md",
             "Usage" => "usage.md",
             "Library" => Any[
                 "Public" => "lib/public.md",
                 "Internals" => "lib/internals.md",
                 "Example matrix generators" => "lib/demos.md"
             ]
         ]
         )
#=
deploydocs(repo = "github.com/RalphAS/Pseudospectra.jl.git",
           osname = "linux",
           julia = "0.6",
           target = "build",
           deps = nothing,
           make = nothing)
=#
