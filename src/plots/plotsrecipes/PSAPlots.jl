module PSAPlots
    using RecipesBase
    using Pseudospectra
    using LinearAlgebra

    import Pseudospectra: vec2ax, _basic_psa_opts, psa_compute

    include("spectralportrait.jl")

    export spectralportrait
end