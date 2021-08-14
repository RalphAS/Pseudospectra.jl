@userplot SpectralPortrait

RecipesBase.@recipe function f(p::SpectralPortrait; npts::Integer = 100)

    if length(p.args) != 1 || !(typeof(p.args[1]) <: AbstractMatrix)
        error("spectralportrait must be given an AbstractMatrix. Got $(typeof(p.args))")
    end
    A₀ = p.args[1]

    local ps_data
    try
        ps_data = new_matrix(A₀)
    catch JE
        @warn "spectralportrait only works for simple cases."
        rethrow(JE)
    end

    A = ps_data.matrix
    ps_dict = ps_data.ps_dict
    B = get(ps_dict, :matrix2, I)
    eigA = ps_dict[:ews]
    zoom = ps_data.zoom_list[ps_data.zoom_pos]
    if isempty(zoom.ax)
        zoom.ax = vec2ax(eigA)
    end
    psa_opts = _basic_psa_opts(zoom, ps_dict)

    Z, xs, ys, t_levels, err, Tproj, eigAproj, algo = psa_compute(A, npts, zoom.ax, eigA, psa_opts, B)

    @series begin
        seriestype := :contour
        xs, ys, log10.(Z)
    end

    @series begin
        seriestype := :scatter
        markercolor := :black
        markersize := 2
        label := "eigvals"
        real.(eigA), imag.(eigA)
    end

end