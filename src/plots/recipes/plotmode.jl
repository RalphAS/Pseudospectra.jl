@userplot PlotMode

function get_mode_title(pseudo, approx, z, s)
    prefix = pseudo ? "Pseudomode: " : "Eigenmode: "
    rel = approx ? "\\approx" : "="
    formatted = @sprintf("%15.2e", s)
    infix = pseudo ? ("\$\\epsilon = " * formatted * "\$") : ("\$\\kappa(\\lambda) " * rel * formatted * "\$")
    λ_str = @sprintf("%12g%+12gi", real(z), imag(z))
    return prefix * infix * "\n  \$\\lambda = " * λ_str * "\$"
end

"""
    plotmode(A, z, U=I; ...)

Plot the (pseudo)eigenmode of `A` associated with (pseudo)eigenvalue `z`.

# Arguments
* `A` an arbitrary square matrix
* `z` the real or complex (pseudo)eigenvalue whose eigenmode is to be plotted.
* `U` optional Schur factor. The eigenmode is premultiplied by `U` before plotting.
# Keyword arguments
* `pseudo::Bool = false` Whether to plot the eigenmode or pseudomode.
* `approx::Bool = false`
* `verbosity::Int = 0` Verbosity parameter passed when computing the eigenmode.

# Examples
```jldoctest
julia> A = Pseudospectra.grcar(40)
julia> plotmode(A, 0.5+2.0im)
```
"""
function plotmode end

RecipesBase.@recipe function f(p::PlotMode, U = I; pseudo=false, approx=false, verbosity=0)
    if length(p.args) != 2
        error("plotmode must be given 2 arguments")
    end
    A = p.args[1]
    z = p.args[2]

    SS, q = get_mode(A, z, pseudo, verbosity)

    q = U * q
    
    #showlog = sum(abs.(diff(q))) * (1/length(q)) >= 10*eps()

    # NB. The title is broken
    #title_str = get_mode_title(pseudo, approx, z, SS[end])
    #title := title_str
    
    @series begin
        seriestype := :path
        linecolor := (pseudo ? :magenta : :cyan)
        label := "Real part"
        real.(q)
    end

    seriestype := :path
    linecolor := :black
    linestyle := :dash

    @series begin
        label := "Modulus"
        abs.(q)
    end

    @series begin
        label := ""
        -abs.(q)
    end
end