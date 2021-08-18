@userplot PlotMode

function get_mode_title(pseudo, approx, z, s)
    prefix = pseudo ? "Pseudomode: " : "Eigenmode: "
    rel = approx ? "\\approx" : "="
    formatted = @sprintf("%15.2e", s)
    infix = pseudo ? ("\$\\epsilon = " * formatted * "\$") : ("\$\\kappa(\\lambda) " * rel * formatted * "\$")
    λ_str = @sprintf("%12g%+12gi", real(z), imag(z))
    return prefix * infix * "\n  \$\\lambda = " * λ_str * "\$"
end

RecipesBase.@recipe f(p::PlotMode; pseudo::Bool, approx::Bool, verbosity)

    z = p.args[1]
    A = p.args[2]
    U = p.args[3]

    SS, q = get_mode(A, z, pseudo, verbosity)

    q = U * q
    
    #showlog = sum(abs.(diff(q))) * (1/length(q)) >= 10*eps()
    title_str = get_mode_title(pseudo, approx, z, SS[end])
    @series begin
        seriestype := :path
        linecolor := (pseudo ? :magenta : :cyan)
        label := "realpt"
        title := title_str
        real.(q)
    end
    @series begin
        seriestype := :path
        linecolor := :black
        label := "abs"
        abs.(q)
    end
    @series begin
        seriestype := :path
        linecolor := :black
        label := ""
        -abs.(q)
    end
end