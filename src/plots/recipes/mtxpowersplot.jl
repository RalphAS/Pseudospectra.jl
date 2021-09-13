@userplot MtxPowersPlot

"""
    mtxpowersplot(A, nmax)

Plot the matrix norm of `A^k` for `k` up to `nmax`

# Examples
```jldoctest
julia> opts = Dict{Symbol,Any}(:npts => 40,:ax => [-110,110,-110,110])
julia> A = diagm(1 => 100.0*ones(4)) + diagm(-2 => 1e-4*ones(3))
5×5 Matrix{Float64}:
 0.0     100.0       0.0       0.0    0.0
 0.0       0.0     100.0       0.0    0.0
 0.0001    0.0       0.0     100.0    0.0
 0.0       0.0001    0.0       0.0  100.0
 0.0       0.0       0.0001    0.0    0.0
 julia> ps_data = new_matrix(A,opts)
 [...]
 julia> driver!(ps_data, opts)
 julia> mtxpowersplot(ps_data, nmax=20)
```
"""
function mtxpowersplot end

RecipesBase.@recipe function f(p::MtxPowersPlot;
    nmax::Integer = 50,
    lbmethod = :none,
    lbdk = 0.25
    )
    if length(p.args) != 1 || !(typeof(p.args[1]) <: PSAStruct)
        error("mtxpowersplot must be given a PSAStruct. Got $(typeof(p.args))")
    end
    ps_data::PSAStruct = p.args[1]
    A = ps_data.input_matrix

    powers, transient = get_mtxpowersnorm(A, nmax)

    alp = maximum(abs.(ps_data.ps_dict[:ews]))

    # NB. We need Plots to use @layout
    #layout := @layout [a; b]
    #link := :x # link the x-axis of both subplots

    # normal plot
    @series begin
        seriestype := :path
        subplot := 1
        label := "Matrix norm"
        powers, transient
    end
    @series begin
        seriestype := :path
        subplot := 1
        label := "Spectral Radius"
        powers, alp.^powers
    end
    # log plot
    #=
    @series begin
        seriestype := :path
        subplot := 2
        label := ""
        yscale := :log10
        powers, transient
    end
    @series begin
        seriestype := :path
        subplot := 2
        label := ""
        yscale := :log10
        powers, alp.^powers
    end
    =#
end