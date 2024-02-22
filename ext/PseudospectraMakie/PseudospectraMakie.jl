#=
 Makie drivers etc. for Pseudospectra.jl

This file is part of Pseudospectra.jl.

Julia implementation
Copyright (c) 2020-2021 Ralph A. Smith

Portions derived from EigTool:
 Copyright (c) 2002-2014, The Chancellor, Masters and Scholars
 of the University of Oxford, and the EigTool Developers. All rights reserved.
 EigTool is maintained on GitHub:  https://github.com/eigtool

SPDX-License-Identifier: BSD-3-Clause
License-Filename: LICENSES/BSD-3-Clause_Eigtool
=#
module PseudospectraMakie

isdefined(Base, :get_extension) ? (using Makie) : (using ..Makie)
using PlotUtils
using Pseudospectra, LinearAlgebra, Printf
# using Colors: RGBA, RGB
using PlotUtils: RGBA, RGB

export MakieGUIState

# we implement specific methods for these here:
import Pseudospectra: redrawcontour, surfplot, arnoldiplotter!, ewsplotter
import Pseudospectra: plotmode, replzdlg, addmark, fillopts, isheadless
import Pseudospectra: mtxexpsplot, mtxpowersplot
import Pseudospectra: zoomin!, zoomout!, _portrait

# we use these internals here:
using Pseudospectra: vec2ax, expandlevels, oneeigcond, psmode_inv_lanczos
using Pseudospectra: dummyqdlg, replqdlg, transient_bestlb
using Pseudospectra: numrange!

const default_opts = Dict{Symbol,Any}(
    :contourstyle => nothing, :markeig => true,
    :no_waitbar => false,
    # No satisfactory colormap choices are built into Plots.
    # So we bit the bullet and implemented the one from EigTool
    # (which is really good); its spec is added below.
    :contourkw => Dict{Symbol,Any}(:linewidth => 3,
                                   :fill => false)
)

struct MakiePlotter <: Pseudospectra.PSAPlotter end

mutable struct MakieGUIState <: GUIState
    mainph # opaque backend object for main plot
    mainfignum::Int
    drawcmd # function to display a plot object (pluggable for GUI use)
    ph2 # opaque backend object for secondary plot
    markerlist
    counter
    is_headless::Bool
    do_savefig::Bool
    figfile

    function MakieGUIState(ph=nothing,num=0,specialcmd=nothing;
                           headless=false, savefigs=true)
        if specialcmd === nothing
            dc = headless ? dcheadless : dcinteractive
        else
            dc = specialcmd
        end
        new(ph,num,dc,nothing,[],0,headless,savefigs,"")
    end
end

isheadless(gs::MakieGUIState) = gs.is_headless

function fillopts(gs::MakieGUIState,optsin::Dict{Symbol,Any}=Dict{Symbol,Any}())
    opts = merge(default_opts,optsin)
end

function drawp(gs::MakieGUIState, p, id)
    gs.drawcmd(gs,p,id)
end

isIJulia() = isdefined(Main, :IJulia) && Main.IJulia.inited

function dcinteractive(gs::MakieGUIState,p,n)
    # CHECKME: do we need special treatment for
    # if isIJulia()
    display(p)
    nothing
end

function dcheadless(gs::MakieGUIState,p,n)
    if gs.do_savefig
        if length(gs.figfile) > 0
            fn = gs.figfile
        else
            fn = tempname() * ".png"
        end
        Makie.save(fn, p)
        @info("graph saved in $(fn)")
    end
    nothing
end

# Makie weirdness
# markersize is in data units (pace Makie docs) so we have to adjust it
const MARKER_RATIO = Ref(0.03)

################################################################
# The main plotting functions

function redrawcontour(gs::MakieGUIState, ps_data::PSAStruct, opts)
    zoom = ps_data.zoom_list[ps_data.zoom_pos]
    ps_dict = ps_data.ps_dict
    levels = expandlevels(zoom.levels)
    Z,x,y = zoom.Z, zoom.x, zoom.y
    if !zoom.computed || isempty(zoom.Z)
        error("should not get here w/ uncomputed zoom")
    end
    eigA = ps_dict[:proj_ews]
    fig = Figure(resolution = (500, 500))
    ax1 = Axis(fig)
    ctx = ax1
    gs.mainph = fig

    if isempty(levels)
        clines = contour!(ctx,x,y,log10.(Z');opts[:contourkw]...)
    else
        kwargs = merge(Dict(:levels => levels),opts[:contourkw])
        clines = contour!(ctx,x,y,log10.(Z'); kwargs...)
        # gs.mainph = quantour(x,y,log10.(Z), levels; opts[:contourkw]...)
    end
    cbar = Colorbar(fig, clines, label="log10(ϵ)", labelpadding=0)
        if !isempty(eigA)
            scatter!(ctx,real(eigA),imag(eigA),color=:black)
        end
        if get(opts,:showimagax,false)
            lines!(ctx,[0,0],zoom.ax[3:4],color=:black, linestyle=:dash)
        end
        if get(opts,:showunitcircle,false)
            lines!(ctx,cos.((0:0.01:2)*π),sin.((0:0.01:2)*π),
                  color=:black, linestyle=:dash)
        end
        if get(opts,:showfov,false)
            if isempty(get(ps_dict,:schur_mtx,[]))
                @warn("showfov set in opts but unavailable")
            else
                if isempty(get(ps_dict,:fov,[]))
                    numrange!(ps_data,get(opts,:fov_npts,20))
                end
                lines!(ctx,real.(ps_dict[:fov]),imag.(ps_dict[:fov]),
                       color=:black)
            end
        end
    setxylims!(ctx,zoom.ax)
    fig[1,1] = ax1
    fig[1,2] = cbar
    drawp(gs,gs.mainph,1)
    nothing
end

"""
make a surface plot of the current spectral portrait
"""
function surfplot(gs::MakieGUIState, ps_data::PSAStruct, opts)
    zoom = ps_data.zoom_list[ps_data.zoom_pos]
    zoom.computed || return
    nx,ny = size(zoom.Z)
    # TODO: allow for user options
    # TODO: control plotting of mesh, e.g.
    # line_freq = floor(Int,min(nx,ny)/10)
    # TODO: add cylinders for eigenvalues, like upstream
    ph = Scene()
    surface!(ph,zoom.x,zoom.y,-log10.(zoom.Z'))
    drawp(gs,ph,2)
end

struct ArnoldiPlotterState{T}
    axis::Scene
    ewsnode::Observable{T}
    shiftsnode::Observable{T}
end

function arnoldiplotter!(gs::MakieGUIState,old_ax,opts,dispvec,infostr,
                         ews,shifts,state)
    update_ax = false

    if !get(opts,:ARPACK_auto_ax,true)
        @warn("non-auto axes for ARPACK not implemented")
        opts[:ARPACK_auto_ax] = true
    end
    if true # get(opts,:ARPACK_auto_ax,true)
        ax = vec2ax(dispvec)
        if  !(state === nothing)
            if isempty(old_ax)
                append!(old_ax,getxylims(state.axis))
                update_ax = true
            else
#                copy!(old_ax,getxylims(state.axis))
                ax[1] = min(ax[1],old_ax[1])
                ax[3] = min(ax[3],old_ax[3])
                ax[2] = max(ax[2],old_ax[2])
                ax[4] = max(ax[4],old_ax[4])
                update_ax = !all(ax .== old_ax)
            end
        else
            update_ax = true
        end
    end # ax setting

    if state === nothing
        ewsnode = Observable(Point2f.(real.(ews),imag.(ews)))
        shiftsnode = Observable(Point2f.(real.(shifts),imag.(shifts)))

        ax1 = Scene()
        ms = MARKER_RATIO[] * (ax[2] - ax[1])
        s1 = scatter!(ax1,ewsnode,color=:black,markersize=ms)
        title(ax1, infostr)

        s2 = scatter!(ax1,shiftsnode,
                      color=:red,marker=:+,markersize=ms,strokewidth=0.1*ms)
        setxylims!(ax1,ax)
        leg = legend([s1,s2],["ews","shifts"])
        scene = Scene()
        vbox(ax1,leg; parent=scene)

        state = ArnoldiPlotterState(ax1, ewsnode, shiftsnode)
        gs.mainph = scene

        # drawcmd() should force a redraw
        drawp(gs,gs.mainph,1)
    else
        state.ewsnode[] = Point2f.(real.(ews),imag.(ews))
        state.shiftsnode[] = Point2f.(real.(shifts),imag.(shifts))
        if update_ax
            ms = MARKER_RATIO[] * (ax[2] - ax[1])
            for j in 1:length(state.axis)
                if state.axis[j] isa Scatter
                    state.axis[j].attributes[:markersize] = ms
                end
            end
            setxylims!(state.axis, ax)
            copy!(old_ax, ax)
        end
    end
    sleep(0.05)
    return state
end

function ewsplotter(gs::MakieGUIState, ews::Vector, zoom)
    if isIJulia()
        # unless we can figure out how to overplot
        return nothing
    end

    fig = Figure()
    gs.mainph = fig
    ax1 = fig[1,1] = Axis(fig)
    scatter!(ax1, real(ews),imag(ews),color=:black)
    setxylims!(ax1,zoom.ax)
    drawp(gs,gs.mainph,1)
end

function plotmode(gs::MakieGUIState,z,A,U,pseudo::Bool,approx::Bool,verbosity)
    m,n = size(A)
    if pseudo
        modestr = "Pseudo"
        niter = 20
        q = normalize!(randn(n) + randn(n)*im)
        Ss,q = psmode_inv_lanczos(A,q,z,1e-15,niter)
        the_col = :magenta
    else
        modestr = "Eigen"
        Ss,q = oneeigcond(A,z,verbosity,
                          dlg=ifelse(isheadless(gs),dummyqdlg,replqdlg))
        the_cond = Ss
        the_col = :cyan
    end
    # FIXME: bail out on error

    # transform to account for Schur decomposition
    q = U * q

    # upstream tries to fetch psmode_x_points_ from base workspace
    x = collect(1:length(q))
    showlog = sum(abs.(diff(q))) * (1/length(q)) >= 10*eps()
    if pseudo
        the_str = @sprintf("ϵ = %15.2e",Ss[end])
    else
        the_rel = approx ? "≈" : "="
        the_str = "κ(λ) " * the_rel *
            @sprintf("%15.2e",the_cond)
    end
    λ_str = @sprintf("%12g%+12gi",real(z),imag(z))
    ax1_title = "$(modestr)mode: " * the_str * "\nλ=" * λ_str
    fig = Figure(resolution = (500, 500))
    ax1 = fig[1,1] = Axis(fig, title = ax1_title)
    l1=lines!(ax1,x,real.(q),color=the_col)
    l2=lines!(ax1,x,abs.(q),color=:black)
    l3=lines!(ax1,x,-abs.(q),color=:black)
    xlims!(ax1,(x[1],x[end]))
    # this should not be needed
    ylims!(ax1,(-1.0,1.0))
    leg = fig[1,2] = Legend(fig, [l1,l2],["realpt", "abs"])


    if showlog
        # TODO: log ticks
        ax2, l2b = lines(fig[2,1],x,log10.(abs.(q)),color=:black)
        xlims!(ax2,x[1],x[end])
        # this should not be needed
        ylims!(ax2,(-10.0,0.0))
        ax2.ylabel = "log₁₀(|xⱼ|)"
        ax2.xlabel = "j"
    else
        ax1.xlabel = "j"
    end
    drawp(gs, fig, 2)
    nothing
end

function mtxpowersplot(gs::MakieGUIState,ps_data::PSAStruct,nmax=50;
                       gradual=false, lbmethod=:none, lbdk=0.25)

    stop_trans = false

    sp_step = 1 # basic stepsize

    A = ps_data.input_matrix^sp_step

    trans = zeros(nmax+1)
    the_pwr = similar(trans)

    trans[1] = 1.0
    the_pwr[1] = 0.0
    ax = [0,20,0,1.5]
    ax2 = copy(ax)
    ax2[3] = 1e-2
    pos = 2
    mtx = Matrix(1.0I,size(A)...)
    max_tp = -Inf
    min_tp = Inf
    ax_factor = 3

    # spectral radius
    alp = maximum(abs.(ps_data.ps_dict[:ews]))
    newfig = true
    pts,bnds = zeros(0),zeros(0)

    function doplot()
        fig = Figure(resolution = (500, 500))
        p1 = fig[1,1] = Axis(fig)
        l1a = lines!(p1,the_pwr[1:pos],trans[1:pos],color=:blue) # dot-dash
        l2a = lines!(p1,the_pwr[1:pos],alp.^the_pwr[1:pos],color=:black, linestyle=:dash)
        xlims!(p1,(ax[1],ax[2]))
        ylims!(p1,(ax[3],ax[4]))
        p1.ylabel = "∥Aᵏ∥"
        p1.xlabel = "k"
        p2, l1b = lines(fig[2,1], the_pwr[1:pos], log10.(trans[1:pos]), color=:blue)
        xlims!(p2,(ax2[1],ax2[2]))
        ylims!(p2,(log10(ax2[3]),log10(ax2[4])))
        p2.ylabel = "log₁₀"
        # TODO: log ticks

        # spectral growth
        l2b = lines!(p2,the_pwr[1:pos],log10.(alp.^the_pwr[1:pos]), color=:black, linestyle=:dash)

        if isempty(bnds)
            leg = fig[3,1] = Legend(fig, [l1b,l2b],["norms", "spectral"])
        else
            l3b = lines!(p2,pts,bnds, linestyle=:dashdot, color=:green)
            leg = fig[3,1] = Legend(fig, [l1b,l2b,l3b],["norms", "spectral", "bound"])
        end
        leg.tellheight = true
        leg.orientation = :horizontal
        # yspan = log10(ax2[4]) - log10(ax2[3])
        # step = max(1,floor(yspan/3))
        #    plot!(yticks = log10(ax2[3]):step:floor(log10(ax2[4])))
        drawp(gs,fig,2)
    end
    if lbmethod != :none
        pts,bnds,sel_pt = transient_bestlb(ps_data,:powers,0:nmax,
                                           method=lbmethod, dk=lbdk)
    end
    while !stop_trans
        mtx = A * mtx
        trans[pos] = norm(mtx)
        the_pwr[pos] = sp_step*(pos-1)
        max_tp = max(max_tp,trans[pos])
        min_tp = min(min_tp,trans[pos])

        if the_pwr[pos] > ax[2]*0.9
            ax[2] *= ax_factor
            ax2[2] *= ax_factor
        end
        if max_tp > ax[4] * 0.9
            ax[4] = max(max_tp * 1.5,ax_factor * ax[4])
        end
        if max_tp > ax2[4] * 0.9
            ax2[4] = ax2[4]^ax_factor
        end
        if (log10(min_tp)-log10(ax2[3])) / (log10(ax2[4]) - log10(ax2[3])) < 0.1
            ax2[3] = ax2[3]^ax_factor
        end

        if gradual
            doplot()
            sleep(0.02)
        end
        if (max_tp > 1e130) || ((min_tp < 1e-130) && (min_tp != 0))
            @warn("stopping: numbers going out of range")
            stop_trans = true
        end
        if min_tp == 0
            @warn("stopping: exactly zero matrix")
            stop_trans = true
        end
        (pos > nmax) && (stop_trans = true)
        pos += 1
    end
    if !gradual
        pos -= 1
        doplot()
    end
end

function mtxexpsplot(gs::MakieGUIState, ps_data::PSAStruct,dt=0.1,nmax=50;
                     gradual=false, lbmethod=:none)

    stop_trans = false
    A = ps_data.input_matrix
    eAdt = exp(dt*A)

    trans = zeros(nmax+1)
    the_time = similar(trans)

    trans[1] = 1.0
    the_time[1] = 0.0
    ax = [0,20*dt,0,2]
    ax2 = copy(ax)
    ax2[3] = 1e-2
    pos = 2
    eAt = Matrix(1.0I,size(A)...)
    max_tp = -Inf
    min_tp = Inf
    ax_factor = 3

    # spectral abscissa
    alp = maximum(real.(ps_data.ps_dict[:ews]))
    newfig = true
    if lbmethod != :none
        pts,bnds,sel_pt = transient_bestlb(ps_data,:exp,(0:nmax)*dt,
                                           method=lbmethod)
    else
        pts,bnds = zeros(0),zeros(0)
    end


    function doplot()
        fig = Figure(resolution = (500, 500))
        p1 = fig[1,1] = Axis(fig)
        l1a = lines!(p1,the_time[1:pos], trans[1:pos],color=:blue)
        l2a = lines!(p1,the_time[1:pos], exp.(alp*the_time[1:pos]), color=:black,
                     linestyle=:dash)
        xlims!(p1,(ax[1],ax[2]))
        ylims!(p1,(ax[3],ax[4]))
        p1.ylabel = "∥exp(At)∥"
        p1.xlabel = "t"

        p2, l1b = lines(fig[2,1], the_time[1:pos], log10.(trans[1:pos]), color=:blue)
        xlims!(p2,(ax2[1],ax2[2]))
        ylims!(p2,(log10(ax2[3]),log10(ax2[4])))
        p2.ylabel = "log₁₀"

        # TODO: logarithmic ticks for y-axis
        l2b = lines!(p2,the_time[1:pos],(1/log(10.0)) * (alp*the_time[1:pos]), color=:black,
                     linestyle=:dash)
        if isempty(bnds)
            leg = fig[3,1] = Legend(fig, [l1b,l2b],["norms", "spectral"])
        else
            l3b = lines!(p2,pts[1:pos],log10.(bnds[1:pos]), color=:green, linestyle=:dashdot)
            leg = fig[3,1] = Legend(fig, [l1b,l2b,l3b],["norms", "spectral", "bound"])
        end
        leg.tellheight = true
        leg.orientation = :horizontal

        # yspan = log10(ax2[4]) - log10(ax2[3])
        #    step = max(1,floor(yspan/3))
        #    plot!(yticks = log10(ax2[3]):step:floor(log10(ax2[4])))
        drawp(gs,fig,2)
    end

    while !stop_trans
        eAt = eAt * eAdt
        trans[pos] = norm(eAt)
        the_time[pos] = (pos-1)*dt
        max_tp = max(max_tp,trans[pos])
        min_tp = min(min_tp,trans[pos])

        if the_time[pos] > ax[2]*0.9
            ax[2] *= ax_factor
            ax2[2] *= ax_factor
        end
        if max_tp > ax[4] * 0.9
            ax[4] = max(max_tp * 1.5,ax_factor * ax[4])
        end
        if max_tp > ax2[4] * 0.9
            ax2[4] = ax2[4]^ax_factor
        end
        if (log10(min_tp)-log10(ax2[3])) / (log10(ax2[4]) - log10(ax2[3])) < 0.1
            ax2[3] = ax2[3]^ax_factor
        end
        if gradual
            doplot()
            sleep(0.02)
        end

        if (max_tp > 1e130) || ((min_tp < 1e-130) && (min_tp != 0))
            @warn("stopping: numbers going out of range")
            stop_trans = true
        end
        (pos > nmax) && (stop_trans = true)
        pos += 1
    end
    if !gradual
        pos -= 1
        doplot()
    end
end

function zoomin!(gs::MakieGUIState, ps_data::PSAStruct,
                 optsin=Dict{Symbol,Any}(); zkw=zeros(0))
    opts = fillopts(gs,optsin)
    if isempty(zkw)
        # User might get z from some other figure, but so what?
        z = replzdlg(gs)
    else
        z = zkw
    end
    ax = getxylims(gs.mainph)
    res = zoomin!(ps_data,z,ax)
    if res < 0
        @warn("unable to zoom based on requested point")
    elseif res == 0
        # redraw existing portrait
        redrawcontour(gs, ps_data, opts)
    else
        # compute and draw new portrait
        driver!(ps_data, opts, gs)
    end
    res
end

function zoomout!(gs::MakieGUIState, ps_data::PSAStruct,
                 optsin=Dict{Symbol,Any}(); include_fov=false, zkw=zeros(0))
    opts = fillopts(gs,optsin)
    if isempty(zkw)
        z = replzdlg(gs)
    else
        z = zkw
    end
    ax = getxylims(gs.mainph); # current plot axes
    res = zoomout!(ps_data,z,ax,include_fov=include_fov)
    if res < 0
        @warn("unable to zoom based on requested point")
    elseif res == 0
        # redraw existing portrait
        redrawcontour(gs, ps_data, opts)
    else
        # compute and draw new portrait
        driver!(ps_data, opts, gs)
    end
    res
end

################################################################
# Wrappers

function _portrait(::MakiePlotter,xs,ys,Z,eigA)
    fig = Figure()
    # ax1 = Axis(fig[1,1])
    # c = contour!(ax1,xs,ys,log10.(Z'))
    ax1, c = contour(fig[1,1],xs,ys,log10.(Z'))
    cbar = Colorbar(fig[1,2], c)
    # ax = getxylims()
    # ms = MARKER_RATIO[] * (ax[2] - ax[1])
    scatter!(ax1, real(eigA), imag(eigA), color=:black, # label="eigvals",
             )
             # markersize=ms)
    # fig[1,1] = ax1
    # fig[1,2] = cbar
    fig
end

"""
    psa(A [,opts]) -> ps_data, graphics_state

Compute and plot pseudospectra of a matrix.

This is a rudimentary command-line driver for the Pseudospectra package.
"""
function psa(::MakiePlotter, A::AbstractMatrix, opts::Dict{Symbol,Any}=Dict{Symbol,Any}())
    allopts = merge(default_opts,optsin)
    ps_data = new_matrix(A,allopts)
        gs = MakieGUIState()

    driver!(ps_data,allopts,gs)
    ps_data, gs
end
#############################################################
# Utilities

function setxylims!(ph,ax)
    xlims!(ph,(ax[1],ax[2]))
    ylims!(ph,(ax[3],ax[4]))
end

function getxylims(ph)
    if ph isa Scene
        sl = scene_limits(ph)
    else
        ctx = content(ph[1,1])
        sl = ctx.finallimits[]
    end
    if sl === nothing
        return fill(0.0,4)
    end
    v0 = sl.origin
    v1 = v0 + sl.widths
    return [v0[1],v1[1],v0[2],v1[2]]
end

# this is here because other backends allow picking with cursor
function replzdlg(gs::MakieGUIState; what="a z-value")
    x,y = NaN,NaN

    println("Specify " * what *
            "; enter real and imaginary components separated by whitespace")
    print("z = ")
    l = readline()
    strs = split(l)
    try
        x = parse(Float64,strs[1])
        y = parse(Float64,strs[2])
    catch
        println("unable to parse input")
        x = NaN
    end
    z = x + y*im
    return z
end

function addmark(gs::MakieGUIState,z,mykey)
    x,y = real(z),imag(z)
    if mykey == :pseudo
        mkey = :circle
        ckey = :magenta
    elseif mykey == :eigen
        mkey = :circle
        ckey = :cyan
    else
        mkey = :x
        ckey = :black
    end
    # FIXME: add marker to list (replacing if any)
    # then trigger a redraw that handles everything
    ctx = gs.mainph[1,1]
    scatter!(ctx,[x],[y],color=ckey,marker=mkey,markersize=7)
    drawp(gs,gs.mainph,1)
    nothing
end

"""
    The colormap from EigTool. This provides clearly
    distinguishable colors for a line-contour plot.
"""
function etcgrad()
    cm = [
        0   1.00000000000000   0.10000000000000
        0   0.92500000000000   0.17500000000000
        0   0.85000000000000   0.25000000000000
        0   0.77500000000000   0.32500000000000
        0   0.70000000000000   0.40000000000000
        0   0.62500000000000   0.47500000000000
        0   0.55000000000000   0.55000000000000
        0   0.47500000000000   0.62500000000000
        0   0.40000000000000   0.70000000000000
        0   0.32500000000000   0.77500000000000
        0   0.25000000000000   0.85000000000000
        0   0.17500000000000   0.92500000000000
        0   0.10000000000000   1.00000000000000
        0.06153846153846   0.09230769230769   1.00000000000000
        0.12307692307692   0.08461538461538   1.00000000000000
        0.18461538461538   0.07692307692308   1.00000000000000
        0.24615384615385   0.06923076923077   1.00000000000000
        0.30769230769231   0.06153846153846   1.00000000000000
        0.36923076923077   0.05384615384615   1.00000000000000
        0.43076923076923   0.04615384615385   1.00000000000000
        0.49230769230769   0.03846153846154   1.00000000000000
        0.55384615384615   0.03076923076923   1.00000000000000
        0.61538461538462   0.02307692307692   1.00000000000000
        0.67692307692308   0.01538461538462   1.00000000000000
        0.73846153846154   0.00769230769231   1.00000000000000
        0.80000000000000                  0   1.00000000000000
        0.80769230769231                  0   0.92307692307692
        0.81538461538462                  0   0.84615384615385
        0.82307692307692                  0   0.76923076923077
        0.83076923076923                  0   0.69230769230769
        0.83846153846154                  0   0.61538461538462
        0.84615384615385                  0   0.53846153846154
        0.85384615384615                  0   0.46153846153846
        0.86153846153846                  0   0.38461538461538
        0.86923076923077                  0   0.30769230769231
        0.87692307692308                  0   0.23076923076923
        0.88461538461538                  0   0.15384615384615
        0.89230769230769                  0   0.07692307692308
        0.90000000000000                  0                  0
        0.86923076923077                  0                  0
        0.83846153846154                  0                  0
        0.80769230769231                  0                  0
        0.77692307692308                  0                  0
        0.74615384615385                  0                  0
        0.71538461538462                  0                  0
        0.68461538461538                  0                  0
        0.65384615384615                  0                  0
        0.62307692307692                  0                  0
        0.59230769230769                  0                  0
        0.56153846153846                  0                  0
        0.53076923076923                  0                  0
        0.50000000000000                  0                  0
        0.54166666666667   0.05000000000000                  0
        0.58333333333333   0.10000000000000                  0
        0.62500000000000   0.15000000000000                  0
        0.66666666666667   0.20000000000000                  0
        0.70833333333333   0.25000000000000                  0
        0.75000000000000   0.30000000000000                  0
        0.79166666666667   0.35000000000000                  0
        0.83333333333333   0.40000000000000                  0
        0.87500000000000   0.45000000000000                  0
        0.91666666666667   0.50000000000000                  0
        0.95833333333333   0.55000000000000                  0
        1.00000000000000   0.60000000000000                  0
    ]
    nc = size(cm,1)
    cc = Vector{RGBA{Float64}}()
    vv = zeros(nc)
    for i in 1:nc
        col = RGBA(RGB(cm[i,:]...))
        push!(cc,col)
        vv[i]=(i-1)/(nc-1)
    end
    if :ContinuousColorGradient in names(PlotUtils, all=true)
        return PlotUtils.ContinuousColorGradient(cc,vv)
    else
        return PlotUtils.ColorGradient(cc,vv)
    end
end

const eigtool_cgrad = etcgrad()

default_opts[:contourkw][:colormap] = eigtool_cgrad;

function __init__()
    Pseudospectra._register_plotter(:Makie, MakieGUIState, MakiePlotter())
end

end # module
