module PSAPlots
    using RecipesBase
    using Pseudospectra
    using LinearAlgebra
    using PlotUtils
    using Printf

    import Pseudospectra: vec2ax, psa_compute
    import Pseudospectra: psmode_inv_lanczos, oneeigcond

    include("utils.jl")
    include("defaults.jl")
    #include("colormap.jl")
    
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