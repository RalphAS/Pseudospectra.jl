@userplot MtxExpsPlot


RecipesBase.@recipe f(p::MtxExpsPlot, dt=0.1, nmax=50, lbmethod=:none)
    if length(p.args) != 1 || !(typeof(p.args[1]) <: PSAStruct)
        error("mtxexpsplot must be given a PSAStruct. Got $(typeof(p.args))")
    end
    ps_data::PSAStruct = p.args[1] # NB hint the type to the compiler
    A = ps_data.input_matrix

    ts, transient = get_mtxexpnorm(A, dt, nmax)
    alp = maximum(real.(ps_data.ps_dict[:ews]))

    layout := @layout [a; b]
    link := :x
    seriestype := :path
    label := ""

    @series begin
        subplot := 1
        ts, transient
    end
    @series begin
        subplot := 1
        ts, exp.(alp*ts)
    end
    @series begin
        subplot := 2
        yscale := :log10
        ts, transient
    end
    @series begin
        subplot := 2
        ts, exp.(alp*ts)
    end
end