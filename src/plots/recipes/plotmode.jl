@userplot PlotMode

struct PlotModeTheme
    linecolor::Symbol
    prefix::String
    
end



RecipesBase.@recipe f(p::PlotMode; pseudo::Bool, approx::Bool, verbosity)

    z = p.args[1]
    A = p.args[2]
    U = p.args[3]

    SS, q = get_mode(A, z, pseudo, verbosity)

    q = U * q
    
    #showlog = sum(abs.(diff(q))) * (1/length(q)) >= 10*eps()

    @series begin
        seriestype := :path
        label := "realpt"
        title := ""
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