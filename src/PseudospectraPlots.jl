#=
 Plots.jl drivers etc. for Pseudospectra.jl

This file is part of Pseudospectra.jl, whose LICENSE file applies.

Julia implementation
Copyright (c) 2017 Ralph A. Smith

Portions derived from EigTool:
 Copyright (c) 2002-2014, The Chancellor, Masters and Scholars
 of the University of Oxford, and the EigTool Developers. All rights reserved.
 EigTool is maintained on GitHub:  https://github.com/eigtool
=#
module PseudospectraPlots

using Plots
using Pseudospectra

export PlotsGUIState, psa, psasimple

# we implement specific methods for these here:
import Pseudospectra: redrawcontour, surfplot, arnoldiplotter!, ewsplotter
import Pseudospectra: plotmode, replzdlg, addmark, fillopts
import Pseudospectra: mtxexpsplot, mtxpowersplot, isheadless

# we use these internals here:
import Pseudospectra: vec2ax, expandlevels, oneeigcond, psmode_inv_lanczos
import Pseudospectra: dummyqdlg, replqdlg, transient_bestlb
import Pseudospectra: numrange!

# for convenience (so one could import just this module)
import Pseudospectra: new_matrix, driver!
export new_matrix, driver!

const default_opts = Dict{Symbol,Any}(
    :contourstyle => nothing, :markeig => true,
    :no_waitbar => false,
    # No satisfactory colormap choices are built into Plots.
    # So we bit the bullet and implemented the one from EigTool
    # (which is really good); its spec is added below.
    :contourkw => Dict{Symbol,Any}(:linewidth => 3,
                                   :fill => false)
)

type PlotsGUIState <: GUIState
    mainph # opaque backend object for main plot
    mainfignum::Int
    drawcmd # function to display a plot object (pluggable for GUI use)
    ph2 # opaque backend object for secondary plot
    markerlist
    counter
    is_headless::Bool
    do_savefig::Bool
    figfile

    function PlotsGUIState(ph=nothing,num=0,specialcmd=nothing;
                           headless=false, savefigs=true)
        if specialcmd == nothing
            dc = headless ? dcheadless : dcinteractive
        else
            dc = specialcmd
        end
        new(ph,num,dc,nothing,[],0,headless,savefigs,"")
    end
end

isheadless(gs::PlotsGUIState) = gs.is_headless

function fillopts(gs::PlotsGUIState,optsin::Dict{Symbol,Any}=Dict{Symbol,Any}())
    opts = merge(default_opts,optsin)
end

function drawp(gs::PlotsGUIState, p, id)
    gs.drawcmd(gs,p,id)
end

function dcinteractive(gs::PlotsGUIState,p,n)
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        return nothing
    end
    Plots.gui(p)
    nothing
end

function dcheadless(gs::PlotsGUIState,p,n)
    if gs.do_savefig
        if length(gs.figfile) > 0
            fn = gs.figfile
        else
            fn = tempname()
        end
        png(p,fn)
        println("graph saved in $(fn).png")
    end
    nothing
end

################################################################
# The main plotting functions

function redrawcontour(gs::PlotsGUIState, ps_data::PSAStruct, opts)
    zoom = ps_data.zoom_list[ps_data.zoom_pos]
    ps_dict = ps_data.ps_dict
    levels = expandlevels(zoom.levels)
    Z,x,y = zoom.Z, zoom.x, zoom.y
    if !zoom.computed || isempty(zoom.Z)
        error("should not get here w/ uncomputed zoom")
    end
    eigA = ps_dict[:proj_ews]
        if isempty(levels)
            gs.mainph = contour(x,y,log10.(Z);opts[:contourkw]...)
        else
            kwargs = merge(Dict(:levels => levels),opts[:contourkw])
            gs.mainph = contour(x,y,log10.(Z); kwargs...)
            # gs.mainph = quantour(x,y,log10.(Z), levels; opts[:contourkw]...)
        end
        setxylims!(gs.mainph,zoom.ax)
        if !isempty(eigA)
            scatter!(gs.mainph,real(eigA),imag(eigA),color="black",label="")
        end
        if get(opts,:showimagax,false)
            plot!(gs.mainph,[0,0],zoom.ax[3:4],color="black",label="")
        end
        if get(opts,:showunitcircle,false)
            plot!(gs.mainph,cos.((0:0.01:2)*π),sin.((0:0.01:2)*π),
                  color="black",label="")
        end
        if get(opts,:showfov,false)
            if isempty(get(ps_dict,:schur_mtx,[]))
                warn("showfov set in opts but unavailable")
            else
                if isempty(get(ps_dict,:fov,[]))
                    numrange!(ps_data,get(opts,:fov_npts,20))
                end
                plot!(gs.mainph,real.(ps_dict[:fov]),imag.(ps_dict[:fov]),
                      color="black",label="")
            end
        end
        drawp(gs,gs.mainph,1)
    nothing
end

"""
make a surface plot of the current spectral portrait
"""
function surfplot(gs::PlotsGUIState, ps_data, opts)
    zoom = ps_data.zoom_list[ps_data.zoom_pos]
    zoom.computed || return
    nx,ny = size(zoom.Z)
    # TODO: allow for user options
    # TODO: control plotting of mesh, e.g.
    # line_freq = floor(Int,min(nx,ny)/10)
    # TODO: add cylinders for eigenvalues, like upstream
    ph = surface(zoom.x,zoom.y,-log10.(zoom.Z))
    drawp(gs,ph,2)
end

function arnoldiplotter!(gs::PlotsGUIState,old_ax,opts,dispvec,infostr,ews,shifts)
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

    gs.mainph = scatter(real(ews),imag(ews),
                        label="ews",
                        color="black", title=infostr)
    setxylims!(gs.mainph,ax)
    scatter!(gs.mainph,real.(shifts),imag.(shifts),
             label="shifts",
             color="red",marker=(5.0,:+,stroke(1)))
    # drawcmd() should force a redraw
    drawp(gs,gs.mainph,1)
    #   sleep(0.01)
end

function ewsplotter(gs::PlotsGUIState, ews::Vector, zoom)
    gs.mainph = scatter(real(ews),imag(ews),color="black",
                            label="")
    setxylims!(gs.mainph,zoom.ax)
    drawp(gs,gs.mainph,1)
end

function plotmode(gs::PlotsGUIState,z,A,U,pseudo::Bool,approx::Bool,verbosity)
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
    p1=plot(x,real.(q),color=the_col,label="realpt",xlims=(x[1],x[end]),
            overwrite=false)
    plot!(p1,x,abs.(q),color=:black,label="abs")
    plot!(p1,x,-abs.(q),color=:black,label="")
    title!(p1,"$(modestr)mode: " * the_str *
           "\n\$\\lambda=" * λ_str * "\$")
    if showlog
        p2=plot(x,abs.(q),color=:black,label="abs",xlims=(x[1],x[end]),
                yscale=:log10,overwrite=false)
        p3 = plot(p1,p2,layout=(2,1),overwrite=false)
        drawp(gs,p3,2)
    else
        drawp(gs,p1,2)
    end
    nothing
end

function mtxpowersplot(gs::PlotsGUIState,ps_data,nmax=50;gradual=false,
                       lbmethod=:none, lbdk=0.25)

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
    pts,bnds = zeros(0),zeros(0)

    function doplot()
            p1 = plot(the_pwr[1:pos],trans[1:pos],label="",
                      overwrite=false) # ".-"
            plot!(p1,the_pwr[1:pos],alp.^the_pwr[1:pos],label="") # "k--"
            plot!(p1,xlims=(ax[1],ax[2]),ylims=(ax[3],ax[4]))

            p2 = plot(the_pwr[1:pos],trans[1:pos],label="",
                      xlims=(ax2[1],ax2[2]),ylims=(ax2[3],ax2[4]),
                      overwrite=false)
            yaxis!(p2,:log10)
            #              yaxis=(:log10,),xlims=(ax2[1],ax2[2]),
            # spectral growth
            plot!(p2,the_pwr[1:pos],alp.^the_pwr[1:pos],label="spectral")

            if !isempty(bnds)
                plot!(p2,pts,bnds,label="lower bound")
            end
            yspan = log10(ax2[4]) - log10(ax2[3])
            step = max(1,floor(yspan/3))
            #    plot!(yticks = log10(ax2[3]):step:floor(log10(ax2[4])))
            p3 = plot(p1,p2,layout=(2,1),overwrite=false)
            drawp(gs,p3,2)
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

function mtxexpsplot(gs::PlotsGUIState,ps_data,dt=0.1,nmax=50; gradual=false)

    stop_trans = false
    A = ps_data.input_matrix
    eAdt = expm(dt*A)

    trans = zeros(nmax+1)
    the_time = similar(trans)

    trans[1] = 1.0
    the_time[1] = 0.0
    ax = [0,20*dt,0,2]
    ax2 = copy(ax)
    ax2[3] = 1e-2
    pos = 2
    eAt = eye(size(A,1))
    max_tp = -Inf
    min_tp = Inf
    ax_factor = 3

    # spectral abscissa
    alp = maximum(real.(ps_data.ps_dict[:ews]))
    newfig = true

    function doplot()
            p1 = plot(the_time[1:pos],trans[1:pos],label="",
                      overwrite=false) # ".-"
            plot!(p1,the_time[1:pos],exp.(alp*the_time[1:pos]),label="") # "k--"
            plot!(p1,xlims=(ax[1],ax[2]),ylims=(ax[3],ax[4]))

            p2 = plot(the_time[1:pos],trans[1:pos],label="",
                      xlims=(ax2[1],ax2[2]),ylims=(ax2[3],ax2[4]),
                      overwrite=false)
            yaxis!(p2,:log10)
#              yaxis=(:log10,),xlims=(ax2[1],ax2[2]),
            plot!(p2,the_time[1:pos],exp.(alp*the_time[1:pos]),label="")

            yspan = log10(ax2[4]) - log10(ax2[3])
            step = max(1,floor(yspan/3))
        #    plot!(yticks = log10(ax2[3]):step:floor(log10(ax2[4])))
            p3 = plot(p1,p2,layout=(2,1),overwrite=false)
            drawp(gs,p3,2)
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

This is a rudimentary command-line driver for the Pseudospectra package.
"""
function psa(A::AbstractMatrix, optsin::Dict{Symbol,Any}=Dict{Symbol,Any}())
    opts = merge(default_opts,optsin)
    ps_data = new_matrix(A,opts)
        gs = PlotsGUIState()

    driver!(ps_data,opts,gs)
    ps_data, gs
end
#############################################################
# Utilities

function setxylims!(ph,ax)
    plot!(ph,xlims=(ax[1],ax[2]),ylims=(ax[3],ax[4]))
end

function getxylims(ph)
    xa=xlims(ph)
    ya=ylims(ph)
    return [xa[1],xa[2],ya[1],ya[2]]
end

# this is here because other backends allow picking with cursor
function replzdlg(gs::PlotsGUIState; what="a z-value")
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

function addmark(gs::PlotsGUIState,z,mykey)
    x,y = real(z),imag(z)
    if mykey == :pseudo
        mkey = (7,:o)
        ckey = :magenta
    elseif mykey == :eigen
        mkey = (7,:o)
        ckey = :cyan
    else
        mkey = (7,:x)
        ckey = :black
    end
    # FIXME: add marker to list (replacing if any)
    # then trigger a redraw that handles everything
    scatter!(gs.mainph,[x],[y],color=ckey,marker=mkey,label="")
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
    cc = Vector{RGBA{Float64}}(0)
    vv = zeros(nc)
    for i in 1:nc
        col = RGBA(RGB(cm[i,:]...))
        push!(cc,col)
        vv[i]=(i-1)/(nc-1)
    end
    ColorGradient(cc,vv)
end

const eigtool_cgrad = etcgrad()

default_opts[:contourkw][:linecolor] = eigtool_cgrad;

end # module
