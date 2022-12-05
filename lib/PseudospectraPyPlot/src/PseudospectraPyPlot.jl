#=
 PyPlot.jl (Matplotlib) drivers etc. for Pseudospectra.jl

This file is part of Pseudospectra.jl.

Julia implementation
Copyright (c) 2017 Ralph A. Smith

Portions derived from EigTool:
 Copyright (c) 2002-2014, The Chancellor, Masters and Scholars
 of the University of Oxford, and the EigTool Developers. All rights reserved.
 EigTool is maintained on GitHub:  https://github.com/eigtool

SPDX-License-Identifier: BSD-3-Clause
License-Filename: LICENSES/BSD-3-Clause_Eigtool
=#
module PseudospectraPyPlot

using PyPlot
using Pseudospectra, LinearAlgebra, Printf

export MPLGUIState, psa

# we implement specific methods for these here:
import Pseudospectra: redrawcontour, surfplot, arnoldiplotter!, ewsplotter
import Pseudospectra: plotmode, replzdlg, addmark, fillopts, isheadless
import Pseudospectra: mtxexpsplot, mtxpowersplot
import Pseudospectra: zoomin!, zoomout!, _portrait

# we use these internals here:
import Pseudospectra: vec2ax, expandlevels, oneeigcond, psmode_inv_lanczos
import Pseudospectra: dummyqdlg, replqdlg, transient_bestlb
import Pseudospectra: numrange!

# for convenience (so one could import just this module)
import Pseudospectra: new_matrix, driver!
export new_matrix, driver!

const default_opts = Dict{Symbol,Any}(
    :contourstyle => nothing, :markeig => true,
    :project_fraction => 0.5, :scale_equal => false,
    :no_waitbar => false,
    # good colormaps: "jet" w/ alpha=0.7, "Dark2" w/ alpha=0.8
    :contourkw => Dict{Symbol,Any}(:linewidths => 3,
                                                :extend => "both",
                                                :cmap => "jet"),
    :contourfkw => Dict{Symbol,Any}(:alpha => 0.7,
                                                 :extend => "both",
                                                 :cmap => "jet"),
    :cbarkw => Dict{Symbol,Any}(:extendrect => true),
    :fillcontour => true
)

mutable struct MPLGUIState <: GUIState
    mainph # opaque backend object for main plot object
    # the fig. numbers are redundant, but more convenient
    mainfignum::Int
    secondaryph # opaque backend object for secondary plot object
    secondaryfignum::Int
    drawcmd
    markerlist
    counter
    is_headless::Bool
    do_savefig::Bool
    figfile

    function MPLGUIState(ph=nothing,num=0,specialcmd=nothing;
                           headless=false, savefigs=true)
        if specialcmd === nothing
            dc = headless ? dcheadless : dcinteractive
        else
            dc = specialcmd
        end
        if headless # CHECKME
            ioff()
        else
            ion()
        end
        new(ph,num,nothing,0,dc,[],0,headless,savefigs,"")
    end
end

isheadless(gs::MPLGUIState) = gs.is_headless

function fillopts(gs::MPLGUIState, optsin::Dict{Symbol,Any}=Dict{Symbol,Any}())
    opts = merge(default_opts,optsin)
end

function drawp(gs::MPLGUIState, p, id)
    gs.drawcmd(gs,p,id)
end

# isIJulia() = isdefined(Main, :IJulia) && Main.IJulia.inited

function dcinteractive(gs::MPLGUIState,p,n)
    # IJulia might need special handling here
    nothing
end

function dcheadless(gs::MPLGUIState,p,n)
    if gs.do_savefig
        # TODO: make sure p is current
        if length(gs.figfile) > 0
            fn = gs.figfile
        else
            fn = tempname()
        end
        savefig(fn * ".png")
        println("graph saved in $(fn).png")
    end
    nothing
end

################################################################
# Plotting routines

function redrawcontour(gs::MPLGUIState, ps_data::PSAStruct, opts)
    zoom = ps_data.zoom_list[ps_data.zoom_pos]
    ps_dict = ps_data.ps_dict
    levels = expandlevels(zoom.levels)
    Z,x,y = zoom.Z, zoom.x, zoom.y
    if !zoom.computed || isempty(zoom.Z)
        error("should not get here w/ uncomputed zoom")
    end
    eigA = ComplexF64.(ps_dict[:proj_ews])
        figure(gs.mainfignum) # in case the user looked elsewhere
        clf()
        !isempty(eigA) && plot(real(eigA),imag(eigA),"k.")
        setxylims!(gs.mainph,zoom.ax)
        if isempty(levels)
            contour(x,y,log10.(Z);opts[:contourkw]...)
            contourf(x,y,log10.(Z);opts[:contourfkw]...)
        else
            kwargs = merge(Dict(:levels => levels),opts[:contourkw])
            contour(x,y,log10.(Z);kwargs...)
            if get(opts,:fillcontour,true)
                kwargs = merge(Dict(:levels => levels),opts[:contourfkw])
                contourf(x,y,log10.(Z);kwargs...)
            end
        end
        if get(opts,:colorbar,true)
            colorbar(;opts[:cbarkw]...)
        end
        if get(opts,:showimagax,false)
            plot([0,0],zoom.ax[3:4],"k-")
        end
        if get(opts,:showunitcircle,false)
            plot(cos.((0:0.01:2)*π),sin.((0:0.01:2)*π),"k-")
        end
        if get(opts,:showfov,false)
            if isempty(get(ps_dict,:schur_mtx,[]))
                @warn("showfov set in opts but unavailable")
            else
                if isempty(get(ps_dict,:fov,[]))
                    numrange!(ps_data,get(opts,:fov_npts,20))
                end
                plot(real.(ps_dict[:fov]),imag.(ps_dict[:fov]),"k--")
            end
        end
    drawp(gs,gs.mainph,1)
    nothing
end

function surfplot(gs::MPLGUIState, ps_data::PSAStruct, opts)
    zoom = ps_data.zoom_list[ps_data.zoom_pos]
    zoom.computed || return
    nx,ny = size(zoom.Z)
    # TODO: allow for user options
    # TODO: control plotting of mesh, e.g.
    # line_freq = floor(Int,min(nx,ny)/10)
    # TODO: add cylinders for eigenvalues, like upstream
        figure()
        surf(zoom.x,zoom.y,-log10.(zoom.Z))
    drawp(gs,gs.mainph,2)
end

function arnoldiplotter!(gs::MPLGUIState,old_ax,opts,dispvec,infostr,ews,shifts,
                         state)
    if gs.mainph === nothing
        ax = vec2ax(dispvec)
    else
        if get(opts,:ARPACK_auto_ax,true)
            ax = vec2ax(dispvec)
            if isempty(old_ax)
                append!(old_ax,getxylims(gs.mainph))
            else
                copy!(old_ax,getxylims(gs.mainph))
                ax[1] = min(ax[1],old_ax[1])
                ax[3] = min(ax[3],old_ax[3])
                ax[2] = max(ax[2],old_ax[2])
                ax[4] = max(ax[4],old_ax[4])
            end
            # ps_data.zoom_list[ps_data.zoom_pos].ax = ax
        else
            @warn("non-auto axes not implemented")
            # CHECKME: are we sure this is ready?
            # ax = ps_data.zoom_list[ps_data.zoom_pos].ax
        end
    end # ax setting

        # figure(mainfigure)
        print(".") # sometimes apparently needed to force a redraw
        clf()
        # if get(opts,:ProgPSA,false)
        #  TODO: draw pseudospectra
        # end
        plot(real.(shifts),imag.(shifts),"r+")
        plot(real.(ews),imag.(ews),"k.")
        axis(ax)
        #text(?,?,infostr,fontsize=12,fontweight=bold)
        title(infostr)
        drawp(gs,gs.mainph,1)
        draw()
    #   sleep(0.01)
    return nothing
end

"""
    selectfig(gs,main::Bool)

start a new figure, or attach to an existing one.
"""
function selectfig(gs::MPLGUIState,main::Bool)
    if main
        if gs.mainfignum > 0
            figure(gs.mainfignum)
        else
            fh = figure()
            gs.mainph = fh
            gs.mainfignum = fh.number
        end
    else
        if gs.secondaryfignum > 0
            figure(gs.secondaryfignum)
        else
            fh = figure()
            gs.secondaryph = fh
            gs.secondaryfignum = fh.number
        end
    end
end

function ewsplotter(gs::MPLGUIState, ews::Vector, zoom)
    selectfig(gs,true)
    clf()
    plot(real(ews),imag(ews),"k.")
    isempty(zoom.ax) && (zoom.ax = vec2ax(ews))
    setxylims!(gs.mainph,zoom.ax)
    drawp(gs,gs.mainph,1)
end

function plotmode(gs::MPLGUIState,z,A,U,pseudo::Bool,approx::Bool,verbosity)
    m,n = size(A)
    if pseudo
        modestr = "Pseudo"
        niter = 20
        q = normalize!(randn(n) + randn(n)*im)
        Ss,q = psmode_inv_lanczos(A,q,z,1e-15,niter)
        the_col = "magenta"
    else
        modestr = "Eigen"
        Ss,q = oneeigcond(A,z,verbosity,
                          dlg=ifelse(isheadless(gs),dummyqdlg,replqdlg))
        the_cond = Ss
        the_col = "cyan"
    end
    # FIXME: bail out on error

    # transform to account for Schur decomposition
    q = U * q

    # upstream tries to fetch psmode_x_points_ from base workspace
    x = collect(1:length(q))
    showlog = sum(abs.(diff(q))) * (1/length(q)) >= 10*eps()
    if pseudo
        the_str = @sprintf("\$\\epsilon = %15.2e\$",Ss[end])
    else
        the_rel = approx ? "\\approx" : "="
        the_str = "\$\\kappa(\\lambda) " * the_rel *
            @sprintf("%15.2e\$",the_cond)
    end
    λ_str = @sprintf("%12g%+12gi",real(z),imag(z))

    selectfig(gs, false)
    clf()
    showlog && subplot(2,1,1)

    plot(x,abs.(q),color="black",label="abs")
    plot(x,-abs.(q),color="black")
    plot(x,real.(q),color=the_col,label="realpt")
    legend()
    xlim(x[1],x[end])
        # DEVNOTE: upstream puts the_str in the plot area (often obscured)
        #=
        y0,y1 = ylim()
        dx = x[end]-x[1]
        dy = y1-y0
        # upstream also has fontweight="bold", but that's hard w/ math mode
        text(x[1]+dx/40,y0+dy/10,the_str,fontsize=12)
        =#
    title("$(modestr)mode: " *
          "\$\\lambda=" * λ_str * "\$\n" * the_str)
    if showlog
        subplot(2,1,2)
        semilogy(x,abs.(q),"k")
        xlim(x[1],x[end])
    end
    drawp(gs,gs.mainph,2)
end

function mtxpowersplot(gs::MPLGUIState,ps_data::PSAStruct,nmax=50;gradual=false,
                       lbmethod=:none,lbdk=0.25)

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

    function doplot()
        # if plotting gradually, control axis limits to avoid jumpiness
        if newfig
            selectfig(gs,false)
            newfig = false
        end
            subplot(2,1,1)
            gradual && cla()
            plot(the_pwr[1:pos],trans[1:pos])
            plot(the_pwr[1:pos],alp.^the_pwr[1:pos],"g--")
            ylabel("\$\\|A^k\\|\$")
            if gradual
                xlim(ax[1],ax[2])
                ylim(ax[3],ax[4])
            end
            subplot(2,1,2)
            gradual && cla()
            semilogy(the_pwr[1:pos],trans[1:pos])
            plot(the_pwr[1:pos],alp.^the_pwr[1:pos],"g--")
            xlabel("k")
            if gradual
                xlim(ax[1],ax[2])
                ylim(ax2[3],ax2[4])
            end
            drawp(gs,gs.mainph,2)
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
    if lbmethod != :none
        pts,bnds,sel_pt = transient_bestlb(ps_data,:powers,0:nmax,
                                           method=lbmethod, dk=lbdk)

        plot(pts,bnds,"m")
        # upstream has something like
        # figure(gs.mainfignum)
        # scatter(sel_pt[:,1],sel_pt[:,2])
    end
end

function mtxexpsplot(gs::MPLGUIState,ps_data::PSAStruct,dt=0.1,nmax=50;
                     gradual=false, lbmethod=:none)

    stop_trans = false
    A = ps_data.input_matrix
    eAdt = exp(dt*A)

    trans = zeros(nmax+1)
    the_time = similar(trans)

    trans[1] = 1.0
    the_time[1] = 0.0
    ax = [0.0,20.0*dt,0,2.0]
    ax2 = copy(ax)
    ax2[3] = 1e-2
    pos = 2
    eAt = Matrix(1.0I,size(A)...)
    max_tp = -Inf
    min_tp = Inf
    ax_factor = 3

    # get spectral abscissa
    alp = maximum(real.(ps_data.ps_dict[:ews]))
    newfig = true

    if lbmethod != :none
        pts,bnds,sel_pt = transient_bestlb(ps_data,:exp,dt*(0:nmax),
                                           method=lbmethod)
    else
        bnds = zeros(0)
    end
    function doplot()
        if newfig
            selectfig(gs,false)
            clf()
            newfig = false
        end
            subplot(2,1,1)
            gradual && cla()
            plot(the_time[1:pos],trans[1:pos])
            plot(the_time[1:pos],exp.(alp*the_time[1:pos]),"g--")
            if gradual
                xlim(ax[1],ax[2])
                ylim(ax[3],ax[4])
            end
            ylabel("\$\\|e^{At}\\|\$")
            subplot(2,1,2)
            gradual && cla()
            semilogy(the_time[1:pos],trans[1:pos])
            plot(the_time[1:pos],exp.(alp*the_time[1:pos]),"g--")
            xlabel("t")
            if !isempty(bnds)
                plot(pts[1:pos],bnds[1:pos],"m")
            end
            if gradual
                xlim(ax[1],ax[2])
                ylim(ax2[3],ax2[4])
            end
        drawp(gs,gs.mainph,2)
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

function zoomin!(gs::MPLGUIState, ps_data::PSAStruct,
                 optsin=Dict{Symbol,Any}(); zkw=zeros(0))
    opts = fillopts(gs,optsin)
    if isempty(zkw)
        # User might get z from some other figure, but so what?
        z = replzdlg(gs)
    else
        z = zkw
    end
    figure(gs.mainfignum)
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

function zoomout!(gs::MPLGUIState, ps_data::PSAStruct,
                 optsin=Dict{Symbol,Any}(); include_fov=false, zkw=zeros(0))
    opts = fillopts(gs,optsin)
    if isempty(zkw)
        z = replzdlg(gs)
    else
        z = zkw
    end
    figure(gs.mainfignum)
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

"""
    psa(A [,opts]) -> ps_data, graphics_state

Compute and plot pseudospectra of a matrix.
This is a rudimentary driver for the Pseudospectra package.
"""
function psa(A::AbstractMatrix, opts::Dict{Symbol,Any}=Dict{Symbol,Any}())
    allopts = merge(default_opts,opts)
    ps_data = new_matrix(A,allopts)
    fh = figure()
    mainfignum = fh.number
    gs = MPLGUIState(fh,mainfignum)
    driver!(ps_data,allopts,gs=gs)
    ps_data, gs
end

function _portrait(xs,ys,Z,eigA)
    p = contour(xs,ys,log10.(Z),cmap="eigtool")
    plot(real(eigA), imag(eigA), "k.", # label="eigvals",
         markersize=5)
    colorbar()
    return p
end
################################################################
# Utilities

function setxylims!(ph,ax)
    # FIXME: this should insure reference to ph
    axis(ax)
end

function getxylims(ph)
    # FIXME: this should insure reference to ph
    collect(axis())
end

function replzdlg(gs::MPLGUIState; what="a z-value")
    x,y = NaN,NaN
    println("click to select $what in main plot frame")
    figure(gs.mainfignum) # in case the user looked elsewhere
    x,y = ginput()[1]
    z = x + y*im
    return z
end

function addmark(gs::MPLGUIState,z,mykey)
    x,y = real(z),imag(z)
    figure(gs.mainfignum)
    #=
    l = gs.markerlist
    # remove any existing mark of this sort
    # WARNING: logic assumes zero or one mark per class
    for i in eachindex(l)
        if l[i][:class] == mykey
            l[i][:handle][1][:remove]()
            deleteat!(l,i)
            break
        end
    end
    =#
    if mykey == :pseudo
        ltype = "mo"
    elseif mykey == :eigen
        ltype = "co"
    else
        ltype = "k*"
    end
    ph = plot([x],[y],ltype)
    # push!(gs.markerlist,Dict(:class => mykey, :handle => ph))
    nothing
end

function et_cmap()
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
    xvals = collect(range(0.0, stop=1.0, length=nc))
    r = [(xvals[i],cm[i,1],cm[i,1]) for i in 1:nc]
    g = [(xvals[i],cm[i,2],cm[i,2]) for i in 1:nc]
    b = [(xvals[i],cm[i,3],cm[i,3]) for i in 1:nc]
    PyPlot.ColorMap("eigtool",r,g,b,
                    Array{Tuple{Float64,Float64,Float64}}(undef, 0),
                    256,1.0)
end

PyPlot.register_cmap("eigtool",et_cmap())

function __init__()
    Pseudospectra._register_plotter(:PyPlot, MPLGUIState)
end

end # module
