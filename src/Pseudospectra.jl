module Pseudospectra
#=
Eigenvalue and Pseudospectrum Analysis for Julia

The Pseudospectra.jl package is a translation of EigTool, but no endorsement
or promotion by the authors of EigTool is implied.

This package is released under a BSD license, as described in the LICENSE file.

Julia code and supplements
Copyright (c) 2017, Ralph Smith

Portions derived from EigTool:

 Copyright (c) 2002-2014, The Chancellor, Masters and Scholars
 of the University of Oxford, and the EigTool Developers. All rights reserved.
 EigTool is maintained on GitHub:  https://github.com/eigtool
=#

using ProgressMeter
# needed for .& and abstract types
using Compat

export new_matrix, driver!
export psa_compute, pspa_2way, pspr_2way
export modeplot, mtxexpsplot, mtxpowersplot, isheadless
export PSAStruct, ArpackOptions, Portrait, GUIState

# Not exported, but may be used by plotting packages:
# vec2ax, expandlevels
# oneeigcond, psmode_inv_lanczos

# Associated plotting packages should provide these, specialized on their
# own GUIState types:
# redrawcontour, surfplot, arnoldiplotter!, ewsplotter, plotmode,
# replzdlg, addmark,

const smallσ = 1e-150

"""
by default, sparse matrices of this size or smaller are converted to full
for pseudospectra computation.
"""
const nmax4autofull = 200
"""
by default, iterative methods are used for computing pseudospectra of dense
matrices of this size or larger.
"""
const nmin4autoiter = 1600

const myname = "PSA"

"""
object to hold state for the GUI used by Pseudospectra
"""
@compat abstract type GUIState end

type LevelDesc
    isunif::Bool
    full_levels
    first
    step
    last
    function LevelDesc(levs)
        if isa(levs,Range)
            levels = collect(levs)
        elseif isa(levs,Vector) && (eltype(levs) <: Real)
            levels = levs
        else
            throw(ArgumentError("level specification must be a range or real vector"))
        end
        isunif = ((length(levels) == 2)
                  || ((length(levels) > 2)
                      && (maximum(abs.(diff(diff(levels)))) < 10*eps())))
        if isunif
            first = levels[1]
            step = levels[2]-levels[1]
            last = levels[end]
            full_levels = nothing
        else
            full_levels = copy(levels)
            first = nothing; step = nothing; last = nothing
        end
        new(isunif,full_levels,first,step,last)
    end
end

type Portrait
    x::Vector
    y::Vector
    Z::Matrix
    npts::Int
    ax::Vector
    levels::LevelDesc
    autolev::Bool
    proj_lev
    dims
    computed::Bool
    scale_equal::Bool
end

import Base: show

function show(io::IO,z::Portrait)
    k = isempty(z.Z) ? "missing" : "present"
    print(io,"npts: $(z.npts), ax: $(z.ax), autolev: $(z.autolev), "
          * "computed: $(z.computed), data $k, proj $(z.proj_lev), "
          * "scale_eq: $(z.scale_equal)")
    println()
    print(io,"levels: $(z.levels)")
end

type ArpackOptions{T}
    # Control of iterative computations:
    nev::Int # nr. of eigenvalues for eigs() to search for
    ncv::Int # max. subspace size for eigs()
    which::Symbol # as for eigs()
    tol::Real # tolerance for eigs()
    maxiter::Int # bound for eigs()
    v0::Vector # initial Ritz vector for eigs()
    have_v0::Bool
    sigma
    #    function ArpackOptions{T}(; nev=6, ncv=0, which=:LM, tol=zero(T),
    function (::Type{ArpackOptions{T}}){T}(; nev=6, ncv=0, which=:LM,
                                           tol=zero(T),maxiter=300,
                                           v0=Vector{T}(0), have_v0=false,
                                           sigma=nothing)
        have_v0 = !isempty(v0)
        if ncv==0
            ncv = max(20,2*nev+1)
        end
        new{T}(nev,ncv,which,tol,maxiter,v0,have_v0,sigma)
    end
end

"""
Wrapper structure for Pseudospectra session data
"""
type PSAStruct
    matrix
    unitary_mtx
    input_matrix
    input_unitary_mtx
    matrix2
    ps_dict::Dict{Symbol,Any}
    zoom_list::Vector{Portrait}
    zoom_pos::Int
    # aside: why is there a :proj_lev entry in ps_dict?
    function PSAStruct(m1,u1,m1i,u1i,m2,dict)
        new(m1,u1,m1i,u1i,m2,dict,Vector{Portrait}(0),0)
    end
end

# Placeholders for plot-specific code implemented elsewhere
function redrawcontour end
function surfplot end
function arnoldiplotter! end
function ewsplotter end
function plotmode end
function replzdlg end
function addmark end
function mtxexpsplot end
function mtxpowersplot end
function fillopts end
function isheadless end

include("compute.jl")
include("utils.jl")
include("xeigs.jl")
include("modes.jl")
include("abscissa.jl")
include("radius.jl")
include("numrange.jl")
include("plotter.jl")

"""
    recalc_levels(σ, ax)

construct a reasonable set of contour levels for displaying `log10(σ)`
where `σ` is a meshed field on axes `ax`.
"""
function recalc_levels(sigmin,ax)
    err = 0
    smin,smax = extrema(sigmin)
    if smax <= smallσ
        levels = Vector{typeof(smin)}(0)
        err = -2
        return levels,err
    end
    if smin <= smallσ
        smin=minimum(sigmin[sigmin >= smallσ])
    end
    max_val = log10(smax)
    min_val = log10(smin)
    num_lines = 8
    scale = min(ax[4]-ax[3],ax[2]-ax[1])
    last_lev = log10(0.03*scale)
    if max_val >= last_lev
        num_lines = ceil(Int,num_lines*(last_lev-min_val)/(max_val-min_val))
        max_val = last_lev
    end
    if max_val < min_val
        max_val = log10(smin) + 0.1*log10(smax/smin)
        num_lines = 3
    end
    max_lines = num_lines
    if num_lines > 0
        stepsize = max(round((max_val-min_val)/num_lines*4)/4,0.25)
        (stepsize >= 0.75) && (stepsize = 1.0)
        if (max_val - min_val)/stepsize < max_lines
            stepsize = max(round((max_val-min_val)/num_lines*10)/10,0.1)
            (stepsize == 0.3) && (stepsize = 0.2)
            (stepsize == 0.4) && (stepsize = 0.5)
            ((stepsize >= 0.6) && (stepsize <= 0.8)) && (stepsize = 0.5)
            (stepsize >= 0.9) && (stepsize = 1.0)
        end
        if (max_val - min_val)/stepsize < ceil(max_lines/2)
            stepsize = (max_val - min_val) / max_lines
        else
            min_val = ceil(min_val/stepsize) * stepsize
            max_val = floor(max_val/stepsize) * stepsize
        end
        if stepsize > 0
            levels = collect(min_val:stepsize:max_val)
        else
            levels = [min_val,max_val]
        end
    end
    ll = length(levels)
    levels = flipdim(levels[end:-max(floor(Int,ll/9),1):1],1)
    (length(levels) == 1) && (levels = levels .* ones(2))
    # upstream has commented-out normality message logic here
    return levels,err
end

"""
    new_matrix(A::AbstractMatrix, opts::Dict{Symbol,Any}=()) -> ps_data

process a matrix into the auxiliary data structure used by Pseudospectra.

# Options
- `:direct::Bool`: force use of a direct algorithm for full matrix?
- `:sparse_direct::Bool`: force use of a direct algorithm for sparse matrix?
- `:real_matrix::Bool`: treat `A` as unitarily equivalent to a real matrix?
- `:verbosity::Int`: obvious
- `:eigA`: eigenvalues of `A`, if already known
- `:proj_lev`: projection level (see `psa_compute`)
- `:npts`: edge length of grid for computing and plotting
- `:arpack_opts::ArpackOptions`: see type description
- `:levels::Vector{Real}`: contour levels
- `:ax::Vector{Real}(4)`: bounding box for computation
- `:scale_equal::Bool`: force equal axis scales for spectral portraits?
"""
function new_matrix(A::AbstractMatrix,
                    opts::Dict{Symbol,Any}=Dict{Symbol,Any}())
    m,n=size(A)
    (m >= n) || throw(ArgumentError(
        "Only square or tall rectangular matrices are supported."))
    (issparse(A) && (m != n)) && throw(ArgumentError(
        "Only square sparse matrices are supported."))
    (any(isnan.(A)) || any(isinf.(A))) && throw(ArgumentError(
        "Input matrix has infinite or invalid entries."))

    # flag for the M x (M-1) Hessenberg form
    AisHess = ((m == (n+1)) && all([x == 0 for x in tril(A,-2)]))

    # user may specify that A is unitarily equivalent to a real matrix
    # even if it is complex
    Aisreal = get(opts,:real_matrix, !(eltype(A) <: Complex))

    verbosity = get(opts,:verbosity,1)

    if haskey(opts,:direct)
        direct = opts[:direct]
        (verbosity > 0) && println("opts specifies direct = $direct")
        convert2full = false
    else
        if issparse(A)
            convert2full = (n <= nmax4autofull) && !get(opts,:sparse_direct,false)
            direct = get(opts,:sparse_direct,false) || convert2full
        else
            direct = (n <= nmin4autoiter)
            if (verbosity > 0) && (n > nmin4autoiter)
                # Might merit a warning since it is most likely not expected
                # so iteration options are probably inappropriate.
                # For now, unspec. iteration options already => a warning.
                println("defaulting to iterative for large dense mtx")
            end
        end
    end

    # TODO: check for correctness of
    # proj_lev, levels, ax, arpack stuff, etc. from opts

    # define as placeholder if not provided
    eigA = get(opts,:eigA,Vector{complex(eltype(A))}(0))

    input_unitary_mtx = get(opts,:unitary_mtx,I)
    proj_lev = get(opts,:proj_lev,Inf)
    npts = get(opts,:npts,setgridsize(n,24,80,!Aisreal))
    Aissquare = (m == n)

    local Tschur, U
    haveschur = false
    if Aissquare && !issparse(A) && direct
        (verbosity > 1) && (m > 100) &&
            println("Attempting initial decomposition...")
        # If square, dense, & direct, we prefer a Schur factorization.
        # Checking for schurfact! method should work,
        # but that's just asking for surprises. This should be robust.
        try
            # For some reason Julia devs think a real Schur decomp
            # should shadow the true (not real!) thing
            if eltype(A) <: Complex
                Tschur,U,eigA  = schur(A)
            else
                Tschur,U,eigA  = schur(A+zero(eltype(A))*im)
            end
            haveschur = true
        catch JE
            isa(JE,MethodError) || rethrow(JE)
        end
        if !haveschur && isempty(eigA)
            try
                eigA = eigvals(A)
            catch JE
                isa(JE,MethodError) || rethrow(JE)
                warn("Failed to compute eigenvalues; proceeding without.")
                println("Exception was $JE")
            end
        end
        (verbosity > 1) && (m > 100) && println("...done.")
    end
    if haveschur
        ps_dict = Dict{Symbol,Any}(:Aisreal => Aisreal,
                                   :isHessenberg => AisHess,
                                   :sparse_direct => false,
                                   :schur_mtx => Tschur,
                                   :schur_unitary_mtx => U,
                                   :projection_on => true,
                                   :proj_lev => proj_lev,
                                   :ews => eigA)
        ps_data = PSAStruct(Tschur, input_unitary_mtx * U,
                            A, input_unitary_mtx, I, ps_dict)

    elseif issparse(A) || AisHess || !direct || Aissquare
        # sparse, Hessenberg, iterative, or needing special handling
        forcedirect = (get(opts,:sparse_direct,false) || Aissquare)

        ps_dict = Dict{Symbol,Any}(:Aisreal => Aisreal,
                                   :isHessenberg => AisHess,
                                   :sparse_direct => forcedirect,
                                   :projection_on => false,
                                   :proj_lev => Inf,
                                   :ews => eigA)
        ps_data = PSAStruct(A, input_unitary_mtx, A, input_unitary_mtx, I,
                            ps_dict)
        if !direct
            if !haskey(opts,:arpack_opts)
                warn("setting default ARPACK options")
                ps_dict[:arpack_opts] = ArpackOptions{eltype(A)}()
            else
                isa(opts[:arpack_opts],ArpackOptions) || throw(
                    ArgumentError("type of opts[:arpack_options] must "
                                  * "be ArpackOptions"))
                ps_dict[:arpack_opts] = opts[:arpack_opts]
            end
        elseif issparse(A) && convert2full
            (verbosity > 0) &&
                println("converting to full for direct computation")
            Tschur,U,eigA  = schur(full(ps_data.matrix))
            ps_dict[:schur_mtx] = Tschur
            ps_dict[:schur_unitary_mtx] = U
            ps_data.matrix = Tschur
            ps_data.unitary_mtx = ps_data.input_unitary_mtx * Tschur
            ps_dict[:sparse_direct] = false
            ps_dict[:projection_on] = true
            ps_dict[:ews] = eigA
            ps_dict[:orig_ews] = copy(eigA)
        end
    else # dense, non-square (but not Hessenberg), and direct
        rfS,rfT = rect_fact(A)
        ps_dict = Dict{Symbol,Any}(:Aisreal => Aisreal,
                                   :isHessenberg => false,
                                   :sparse_direct => false,
                                   :projection_on => false,
                                   :proj_lev => Inf,
                                   :ews => eigA)
        ps_data = PSAStruct(rfS,input_unitary_mtx,A,input_unitary_mtx,rfT,
                            ps_dict)
    end
    ps_dict[:orig_ews] = eigA
    ps_dict[:ew_estimates] = false

    if !ps_dict[:isHessenberg] && !isa(ps_data.unitary_mtx,UniformScaling)
        if size(ps_data.unitary_mtx,2) ∉ [m,1]
            ps_data.unitary_mtx = I
        end
    end

    lev = get(opts,:levels,zeros(0))
    zoom = Portrait(zeros(0),zeros(0),zeros(0,0),
                    npts, get(opts,:ax,zeros(0)),
                    LevelDesc(lev),
                    isempty(lev), proj_lev,
                    size(ps_data.matrix), false,
                    get(opts,:scale_equal,false))
    push!(ps_data.zoom_list,zoom)
    ps_data.zoom_pos = 1
    # save for use when returning to initial plot
    ps_dict[:init_opts] = deepcopy(zoom)
    ps_dict[:init_direct] = direct
    ps_dict[:direct] = direct
    ps_dict[:verbosity] = verbosity

    # DEVNOTE: if direct, upstream constructs axes and calls origplot/redraw
    return ps_data
end

"""
    new_matrix(A, opts::Dict{Symbol,Any}=()) -> ps_data

process a linear operator object into the auxiliary data structure used by Pseudospectra.

There must be methods with `A` for `eltype`, `size`, and `A_mul_B!`.
"""
function new_matrix(A, opts::Dict{Symbol,Any}=Dict{Symbol,Any}())
    m,n=size(A)
    (m == n) || throw(ArgumentError(
        "Only square linear operators are supported."))
    Aisreal = get(opts,:real_matrix, !(eltype(A) <: Complex))

    verbosity = get(opts,:verbosity,1)
    direct = false
    eigA = get(opts,:eigA,Vector{complex(eltype(A))}(0))

    input_unitary_mtx = get(opts,:unitary_mtx,I)
    npts = get(opts,:npts,setgridsize(n,24,80,!Aisreal))
    proj_lev = get(opts,:proj_lev,Inf)

    ps_dict = Dict{Symbol,Any}(:Aisreal => Aisreal,
                                   :isHessenberg => false,
                                   :sparse_direct => false,
                                   :projection_on => false,
                                   :proj_lev => Inf,
                                   :ews => eigA)
    ps_data = PSAStruct(A, input_unitary_mtx, A, input_unitary_mtx, I,
                        ps_dict)
    if !haskey(opts,:arpack_opts)
        warn("setting default ARPACK options")
        ps_dict[:arpack_opts] = ArpackOptions{eltype(A)}()
    else
        isa(opts[:arpack_opts],ArpackOptions) || throw(
            ArgumentError("type of opts[:arpack_options] must "
                          * "be ArpackOptions"))
        ps_dict[:arpack_opts] = opts[:arpack_opts]
    end
    lev = get(opts,:levels,zeros(0))
    zoom = Portrait(zeros(0),zeros(0),zeros(0,0),
                    npts, get(opts,:ax,zeros(0)),
                    LevelDesc(lev),
                    isempty(lev), proj_lev,
                    size(ps_data.matrix), false,
                    get(opts,:scale_equal,false))
    push!(ps_data.zoom_list,zoom)

    ps_dict[:orig_ews] = eigA
    ps_dict[:ew_estimates] = false
    ps_data.zoom_pos = 1
    # save for use when returning to initial plot
    ps_dict[:init_opts] = deepcopy(zoom)
    ps_dict[:init_direct] = direct
    ps_dict[:direct] = direct
    ps_dict[:verbosity] = verbosity

    return ps_data
end

"""
    driver!(ps_data, opts, gs)

Compute pseudospectra and plot a spectral portrait.

If using an iterative method to get eigenvalues, invokes that first.

# Arguments
- `ps_data::PSAStruct`: ingested matrix, as processed by `new_matrix`
- `gs::GUIState`: object handling graphical output
- `opts::Dict{Symbol,Any}`: options passed to `redrawcontour`, `arnoldiplotter!`
Other options are already stored in `ps_data`.
"""
function driver!(ps_data::PSAStruct, optsin::Dict{Symbol,Any},
                 gs::GUIState; myprintln=println, mywarn=warn)
    # DEVNOTE: mostly corresponds to switch_redraw.m in EigTool
    opts = fillopts(gs,optsin)
    ps_dict = ps_data.ps_dict
    verbosity = get(ps_dict,:verbosity,1)
    if ps_dict[:direct] || get(ps_dict,:proj_valid,false)
        n,m = size(ps_data.matrix)
        A = ps_data.matrix
        B = ps_data.matrix2
        eigA = ps_dict[:ews]
        zoom = ps_data.zoom_list[ps_data.zoom_pos]

        if !isempty(eigA)
            # This sets the default domain for the typical case
            isempty(zoom.ax) && (zoom.ax = vec2ax(eigA))
            if !isheadless(gs)
                # show eigenvalues while waiting
                ewsplotter(gs, eigA, zoom)
            end
        end

        psa_opts = Dict{Symbol,Any}(:ax=>zoom.ax,:npts=>zoom.npts,
                                    :levels=>expandlevels(zoom.levels),
                                    :recompute_levels=>zoom.autolev,
                                    :proj_lev=>zoom.proj_lev,
                                    :scale_equal=>zoom.scale_equal,
                                    :real_matrix=>ps_dict[:Aisreal],
                                    :verbosity=>verbosity)
        ss = size(A)
        Z,x,y,t_levels,err,Tproj,eigAproj,algo = psa_compute(A,eigA,psa_opts,
                                                        myprintln=myprintln,
                                                        mywarn=mywarn)

        # FIXME: handle projection properly
        ps_dict[:proj_ews] = eigAproj

        if err != 0
            warn("PSA computation failed")
            # FIXME: reset GUI if any
            return nothing
        end

        (verbosity > 1) && println("algorithm: ",algo)
        zoom = ps_data.zoom_list[ps_data.zoom_pos]
        if zoom.autolev
            (verbosity > 1) && println("levels: $t_levels")
            zoom.levels = LevelDesc(t_levels)
        end
        zoom.x = x
        zoom.y = y
        zoom.Z = Z
        zoom.computed = true
        zoom.dims = size(ps_data.matrix)

        redrawcontour(gs, ps_data, opts)
    else
        # iterative method
        ps_data.matrix = ps_data.input_matrix
        m,n = size(ps_data.matrix)
        ao = ps_dict[:arpack_opts]

        function xeigsproducer(chnl)
            local ews,H,V
            local nconv,niter,nmult
            try
                ews,v,nconv,niter,nmult,resid,H,V = xeigs(ps_data.matrix,I,chnl,
                                                      nev=ao.nev,ncv=ao.ncv,
                                                      which=ao.which,tol=ao.tol,
                                                      maxiter=ao.maxiter,
                                                      v0=ao.v0,
                                                      sigma=ao.sigma,
                                                      options=opts)
            catch JE
            # FIXME: post a dialog, reset GUI if any
                warn("eigs failed")
                display(JE)
                println()
                stuff = (:failure,nothing)
                if VERSION < v"0.6-"
                    produce(stuff)
                else
                    put!(chnl,stuff)
                end
                return nothing
            end
        end

        local ews,H,V
        local nconv,niter,nmult
        old_ax = zeros(0)
        if VERSION >= v"0.6-"
            chnl = Channel(xeigsproducer)
            xeigs_result = take!(chnl)
        else
            chnl = Task(()->xeigsproducer(true))
            xeigs_result = consume(chnl)
        end
        while xeigs_result[1] ∉ [:finale,:failure]
            the_key,dispvec,the_str,the_ews,the_shifts = xeigs_result
            if !isheadless(gs)
                arnoldiplotter!(gs,old_ax,opts,dispvec,the_str,the_ews,the_shifts)
            end # if gs
            if VERSION >= v"0.6-"
                xeigs_result = take!(chnl)
            else
                xeigs_result = consume(chnl)
            end
        end
        if xeigs_result[1] == :failure
            return nothing
        end
        the_key,ews,v,nconv,niter,nmult,resid,H,V = xeigs_result
        if verbosity > 0
            println("xeigs: $nconv of $(ao.nev) converged in $niter ($nmult MxV)")
        end
        if verbosity > 1
            println("xeigs ews:")
            display(ews); println()
        end
        ps_dict[:ew_estimates] = true
        ps_dict[:proj_matrix] = H
        ps_data.matrix = H
        ps_dict[:isHessenberg] = true
        ps_dict[:proj_unitary_mtx] = ps_data.input_unitary_mtx * V
        ps_data.unitary_mtx = ps_dict[:proj_unitary_mtx]
        ps_dict[:proj_valid] = true
        ps_dict[:ews] = ews
        # reset zoom list
        ps_data.zoom_pos = 1
        deleteat!(ps_data.zoom_list,2:length(ps_data.zoom_list))
        # CHECKME: do we need remove()?
        ps_dict[:mode_markers] = []
        zoom = ps_data.zoom_list[1]
        origax = ps_dict[:init_opts].ax
        if !isempty(origax) # && !ps_dict[:init_direct]
            copy!(zoom.ax,origax)
            # FIXME: needs correct validity check
#        elseif gs.mainph != nothing
#            println("using eigvals for axis limits")
#            copy!(zoom.ax,getxylims(gs.mainph))
        end
        zoom.autolev = ps_dict[:init_opts].autolev
        zoom.levels = deepcopy(ps_dict[:init_opts].levels)

        delete!(ps_dict,:proj_axes)
        delete!(ps_dict,:comp_proj_lev)
        origplot!(ps_data,opts,gs) # WARNING: reentrant
        # TODO: reset GUI
    end
    nothing
end

"""
called for a new matrix. Makes sure zoom list is ok, then redraws (unless
`ax_only`).
"""
function origplot!(ps_data::PSAStruct, opts, gs; ax_only = false)
    ps_data.zoom_pos = 1
    ps_dict = ps_data.ps_dict
    deleteat!(ps_data.zoom_list,2:length(ps_data.zoom_list))
    zoom = ps_data.zoom_list[1]
    if isempty(zoom.ax)
        if isempty(get(ps_dict,:ews,[]))
            warn("origplot called w/o preset axes or eigenvalues")
        else
            zoom.ax = vec2ax(ps_dict[:ews])
        end
    end
    ax_only || driver!(ps_data,opts,gs)
    nothing
end


################################################################
# FIXME: until we think of a better way to handle this:
include("../examples/demo_mtx.jl")

end # module
