#=
 PyPlot.jl (Matplotlib) drivers etc. for Pseudospectra.jl

This file is part of Pseudospectra.jl, whose LICENSE file applies.

Julia implementation
Copyright (c) 2017 Ralph A. Smith

Portions derived from EigTool:
 Copyright (c) 2002-2014, The Chancellor, Masters and Scholars
 of the University of Oxford, and the EigTool Developers. All rights reserved.
 EigTool is maintained on GitHub:  https://github.com/eigtool
=#
module PseudospectraMPL

using PyPlot
using Pseudospectra

export MPLGUIState, psa, psasimple

# we implement specific methods for these here:
import Pseudospectra: redrawcontour, surfplot, arnoldiplotter!, ewsplotter
import Pseudospectra: plotmode, replzdlg, addmark, fillopts, isheadless
import Pseudospectra: mtxexpsplot, mtxpowersplot

# we use these internals here:
import Pseudospectra: vec2ax, expandlevels, oneeigcond, psmode_inv_lanczos
import Pseudospectra: dummyqdlg, replqdlg
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

type MPLGUIState <: GUIState
    mainph # opaque backend object for main plot object
    mainfignum::Int
    drawcmd
    ph2 # opaque backend object for secondary plot object
    markerlist
    counter
    is_headless::Bool
    do_savefig::Bool
    figfile

    function MPLGUIState(ph=nothing,num=0,specialcmd=nothing;
                           headless=false, savefigs=true)
        if specialcmd == nothing
            dc = headless ? dcheadless : dcinteractive
        else
            dc = specialcmd
        end
        if headless # CHECKME
            ioff()
        else
            ion()
        end
        new(ph,num,dc,nothing,[],0,headless,savefigs,"")
    end
end

isheadless(gs::MPLGUIState) = gs.is_headless

function fillopts(gs::MPLGUIState, optsin::Dict{Symbol,Any}=Dict{Symbol,Any}())
    opts = merge(default_opts,optsin)
end

function drawp(gs::MPLGUIState, p, id)
    gs.drawcmd(gs,p,id)
end

function dcinteractive(gs::MPLGUIState,p,n)
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        return nothing
    end
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
    eigA = ps_dict[:proj_ews]
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
                warn("showfov set in opts but unavailable")
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

function surfplot(gs::MPLGUIState, ps_data, opts)
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

function arnoldiplotter!(gs::MPLGUIState,old_ax,opts,dispvec,infostr,ews,shifts)
    if gs.mainph == nothing
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
            warn("non-auto axes not implemented")
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
end

function ewsplotter(gs::MPLGUIState, ews::Vector, zoom)

        figure(gs.mainfignum)
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

        # FIXME: manage figure no
        fh = figure()
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

function mtxpowersplot(gs::MPLGUIState,ps_data,nmax=50;gradual=false)

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
    mtx = eye(size(A,1))
    max_tp = -Inf
    min_tp = Inf
    ax_factor = 3

    # spectral radius
    alp = maximum(abs.(ps_data.ps_dict[:ews]))
    newfig = true

    function doplot()
        # if plotting gradually, control axis limits to avoid jumpiness
            if newfig
                fh = figure()
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
            warn("stopping: numbers going out of range")
            stop_trans = true
        end
        if min_tp == 0
            warn("stopping: exactly zero matrix")
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

function mtxexpsplot(gs::MPLGUIState,ps_data,dt=0.1,nmax=50; gradual=false)

    stop_trans = false
    A = ps_data.input_matrix
    eAdt = expm(dt*A)

    trans = zeros(nmax+1)
    the_time = similar(trans)

    trans[1] = 1.0
    the_time[1] = 0.0
    ax = [0.0,20.0*dt,0,2.0]
    ax2 = copy(ax)
    ax2[3] = 1e-2
    pos = 2
    eAt = eye(size(A,1))
    max_tp = -Inf
    min_tp = Inf
    ax_factor = 3

    # get spectral abscissa
    alp = maximum(real.(ps_data.ps_dict[:ews]))
    newfig = true

    function doplot()
            if newfig
                fh = figure()
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
            warn("stopping: numbers going out of range")
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

################################################################
# Wrappers

"""
    psa(A [,opts]) -> ps_data, graphics_state

Compute and plot pseudospectra of a matrix.
This is the main command-line driver for the Pseudospectra package.
"""
function psa(A::AbstractMatrix, optsin::Dict{Symbol,Any}=Dict{Symbol,Any}())
    opts = merge(default_opts,optsin)
    ps_data = new_matrix(A,opts)
        fh = figure()
        mainfignum = fh[:number]
        gs = MPLGUIState(fh,mainfignum)
    driver!(ps_data,opts,gs,redrawcontour)
    ps_data, gs
end

"""
    psasimple(A,[,opts])

Compute and plot pseudospectra of a dense square matrix.
This is a simple interface w/o all the data structures.
"""
function psasimple(A::AbstractMatrix, opts::Dict{Symbol,Any}=Dict{Symbol,Any}())
    # from new_matrix()
    n,m=size(A)
    @assert n==m "matrix must be square"
    opts[:real_matrix] = !(eltype(A) <: Complex)
    Tschur,U,eigA  = schur(A)

        fh = figure()
        mainfignum = fh[:number]
        gs = MPLGUIState(fh,mainfignum)
        plot(real(eigA),imag(eigA),"k.")

    if haskey(opts,:ax)
        setxylims!(gs.mainph,opts[:ax])
    else
        println("using eigvals for axis limits")
        opts[:ax] = getxylims(gs.mainph)
    end

    Z,x,y,levels,err,Tproj,eigAproj = psa_compute(Tschur,eigA,opts)

    # from redrawcontour()
    figure(gs.mainfignum) # in case the user looked elsewhere
    clf()
    plot(real(eigA),imag(eigA),"k.")
    setxylims!(gs.mainph,opts[:ax])
    if isempty(levels)
        contour(x,y,log10.(Z),extend="both",linewidths=2)
    else
        contour(x,y,log10.(Z),levels=levels,extend="both",linewidths=2)
    end
    if get(opts,:colorbar,true)
        colorbar(extendrect=true)
    end
    nothing
end

################################################################
# Utilities

function setxylims!(ph,ax)
    axis(ax)
end

function getxylims(ph)
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

end # module
