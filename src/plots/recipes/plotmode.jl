@userplot PlotMode

function get_mode_title(pseudo, approx, z, s)
    prefix = pseudo ? "Pseudomode: " : "Eigenmode: "
    rel = approx ? "\\approx" : "="
    formatted = @sprintf("%15.2e", s)
    infix = pseudo ? ("\$\\epsilon = " * formatted * "\$") : ("\$\\kappa(\\lambda) " * rel * formatted * "\$")
    λ_str = @sprintf("%12g%+12gi", real(z), imag(z))
    return prefix * infix * "\n  \$\\lambda = " * λ_str * "\$"
end

RecipesBase.@recipe function f(p::PlotMode; pseudo=false, approx=false, verbosity=0)
    if length(p.args) != 3
        error("plotmode must be given 3 arguments")
    end
    z = p.args[1]
    A = p.args[2]
    U = p.args[3]

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