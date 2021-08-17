module PSAPlots
    using RecipesBase
    using Pseudospectra
    using LinearAlgebra

    import Pseudospectra: vec2ax, _basic_psa_opts, psa_compute

    include("utils.jl")
    include("colormap.jl")
    
    include("recipes/spectralportrait.jl")
    include("recipes/plotmode.jl")
    include("recipes/mtxpowersplot.jl")
    include("recipes/mtxexpsplot.jl")
    include("recipes/surfplot.jl")

    export spectralportrait, spectralportrait!
    export plotmode, plotmode!
    export mtxpowersplot, mtxpowersplot!
    export mtxexpsplot, mtxexpsplot!
    export surfplot, surfplot!
end