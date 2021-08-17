@userplot MtxPowersPlot


RecipesBase.@recipe function f(p::MtxPowersPlot, nmax = 50,
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

    layout := @layout [a; b]
    link := :x # link the x-axis of both subplots

    # normal plot
    @series begin
        seriestype := :path
        subplot := 1
        label := ""
        powers, transient
    end
    @series begin
        seriestype := :path
        subplot := 1
        label := ""
        powers, alp.^powers
    end
    # log plot
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
end