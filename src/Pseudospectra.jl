module Pseudospectra
#=
Eigenvalue and Pseudospectrum Analysis for Julia

The Pseudospectra.jl package is a translation of EigTool, but no endorsement
or promotion by the authors of EigTool is implied.

This package is released under a BSD license, as described in the LICENSE file.

Julia code and supplements
Copyright (c) 2017-2019, Ralph Smith

Portions derived from EigTool:

 Copyright (c) 2002-2014, The Chancellor, Masters and Scholars
 of the University of Oxford, and the EigTool Developers. All rights reserved.
 EigTool is maintained on GitHub:  https://github.com/eigtool

SPDX-License-Identifier: BSD-3-Clause
License-Filename: LICENSES/BSD-3-Clause_Eigtool
=#

using LinearAlgebra, SparseArrays, Arpack, Printf

using ProgressMeter

export new_matrix, driver!, spectralportrait
export psa_compute, psa_radius, psa_abscissa
export numerical_range, numerical_abscissa
export modeplot, mtxexpsplot, mtxpowersplot, isheadless, iscomputed
export PSAStruct, ArpackOptions, Portrait, GUIState

# Plotting packages should probably extend these:
export zoomin!, zoomout!

# Not exported, but may be used by plotting packages:
# vec2ax, expandlevels, isvalidax
# oneeigcond, psmode_inv_lanczos, transient_bestlb, set_method!

# Associated plotting packages should provide these, specialized on their
# own GUIState types:
# redrawcontour, surfplot, arnoldiplotter!, ewsplotter, plotmode,
# replzdlg, addmark

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

include("types.jl")

# Placeholders for plot-specific code implemented elsewhere
function redrawcontour end
function surfplot end
function arnoldiplotter! end
function _portrait end

"""
    ewsplotter(gs::GUIState, ews::Vector, zoom)

plot eigenvalues

So we have something to look at while waiting for the compute engines.
"""
function ewsplotter end
function plotmode end
function replzdlg end

function addmark end
"""
    mtxexpsplot(ps_data,dt=0.1,nmax=50; gs::GUIState = defaultgs(), gradual=false)

plot the evolution of `∥e^(tA)∥`.

This is useful for analyzing linear initial value problems `∂x/∂t = Ax`.
"""
function mtxexpsplot(ps_data::PSAStruct, dt=0.1, nmax=50;
                     gs::GUIState=defaultgs(), kws...)
    mtxexpsplot(gs, ps_data, dt, nmax; kws...)
end

"""
    mtxpowersplot(ps_data, nmax=50; gs::GUIState = defaultgs(), gradual=false)

plot norms of powers of a matrix `∥A^k∥`

This is useful for analyzing iterative linear algebra methods.
"""
function mtxpowersplot(ps_data::PSAStruct, nmax=50;
                       gs::GUIState = defaultgs(), kws...)
    mtxpowersplot(gs, ps_data, nmax=nmax; kws...)
end

function fillopts end
function isheadless end

include("utils.jl")
include("compute.jl")
include("xeigs.jl")
include("modes.jl")
include("abscissa.jl")
include("radius.jl")
include("numrange.jl")
include("transients.jl")
include("plotter.jl")
include("zooming.jl")

"""
    new_matrix(A::AbstractMatrix, opts::Dict{Symbol,Any}=()) -> ps_data

process a matrix into the auxiliary data structure used by Pseudospectra.

# Options
- `:direct::Bool`: force use of a direct algorithm?
- `:keep_sparse::Bool`: use sparse matrix code even if `A` is not large?
- `:real_matrix::Bool`: treat `A` as unitarily equivalent to a real matrix?
- `:verbosity::Int`: obvious
- `:eigA`: eigenvalues of `A`, if already known
- `:proj_lev`: projection level (see `psa_compute`)
- `:npts`: edge length of grid for computing and plotting pseudospectra
- `:arpack_opts::ArpackOptions`: (see type description)
- `:levels::Vector{Real}`: contour levels (if automatic choice is not wanted)
- `:ax::Vector{Real}(4)`: bounding box for computation `[xmin,xmax,ymin,ymax]`
- `:scale_equal::Bool`: force isotropic axes for spectral portraits?
- `:threaded::Bool`: distribute `Z` values over Julia threads?
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
    # Presumably intended for the case where projection is done
    # by a Krylov scheme outside this package.
    AisHess = ((m == (n+1)) && all([x == 0 for x in tril(A,-2)]))

    # User may specify that A is unitarily equivalent to a real matrix
    # even if it is complex
    Aisreal = get(opts,:real_matrix, !(eltype(A) <: Complex))

    verbosity = get(opts,:verbosity,1)

    convert2full = issparse(A) & (n <= nmax4autofull) &
        !get(opts,:keep_sparse,false)
    if haskey(opts,:direct)
        direct = opts[:direct]
    else
        if issparse(A)
            direct = convert2full
        else
            direct = (n <= nmin4autoiter)
            if (verbosity > 0) && (n > nmin4autoiter)
                # Might merit a warning here since it is most likely not
                # expected so iteration options are probably inappropriate.
                # For now, unspec. iteration options already => a warning.
                println("defaulting to iterative for large dense mtx")
            end
        end
    end

    # TODO: check for correctness of
    # proj_lev, levels, ax, arpack stuff, etc. from opts

    # define as placeholder if not provided
    eigA = get(opts,:eigA,Vector{complex(eltype(A))}())

    input_unitary_mtx = get(opts,:unitary_mtx,I)
    proj_lev = get(opts,:proj_lev,Inf)
    npts = get(opts,:npts,setgridsize(n,24,80,!Aisreal))
    Aissquare = (m == n)

    local Tschur, U
    haveschur = false
    if Aissquare && !issparse(A) && direct
        # if small, delay is negligible
        (verbosity > 1) && (m > 100) &&
            println("Attempting initial decomposition...")
        # If square, dense, & direct, we prefer a Schur factorization.
        # Checking for schurfact! method should work,
        # but that's just asking for surprises. This should be robust.
        try
            if eltype(A) <: Complex
                F = schur(A)
            else
                # For some reason Julia devs think a real Schur decomp
                # should shadow the true (not real!) thing
                F = schur(A .+ zero(eltype(A))*im)
            end
            Tschur,U,eigA  = F.T,F.Z,F.values
            haveschur = true
        catch JE
            isa(JE,MethodError) || rethrow(JE)
        end
        if !haveschur && isempty(eigA)
            try
                eigA = eigvals(A)
            catch JE
                isa(JE,MethodError) || rethrow(JE)
                # ME: maybe trap algorithmic errors too; what are they?

                # Warning is needed here since it explains why we need axes
                # for the driver.
                @warn("Failed to compute eigenvalues; proceeding without.")
                if verbosity > 0
                    # If we display(JE) we get the whole damn matrix too
                    println("Exception was method error: ",JE.f,
                            " for ",typeof(A))
                end
            end
        end
        (verbosity > 1) && (m > 100) && println("...done.")
    end
    if haveschur
        ps_dict = Dict{Symbol,Any}(:Aisreal => Aisreal,
                                   :isHessenberg => AisHess,
                                   :schur_mtx => Tschur,
                                   :schur_unitary_mtx => U,
                                   :direct => true,
                                   :projection_on => true,
                                   :proj_lev => proj_lev,
                                   :ews => eigA)
        ps_data = PSAStruct(UpperTriangular(Tschur), input_unitary_mtx * U,
                            A, input_unitary_mtx, ps_dict)

    elseif issparse(A) || AisHess || !direct || Aissquare
        # sparse, Hessenberg, iterative, or needing special handling
        # CHECKME: previously seemed to force
        #   direct |= Aissquare
        #
        ps_dict = Dict{Symbol,Any}(:Aisreal => Aisreal,
                                   :isHessenberg => AisHess,
                                   :projection_on => false,
                                   :proj_lev => Inf,
                                   :ews => eigA)
        ps_data = PSAStruct(A, input_unitary_mtx, A, input_unitary_mtx,
                            ps_dict)
        if !direct
            if !haskey(opts,:arpack_opts)
                @info("setting default ARPACK options")
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
            Atmp = full(ps_data.matrix)
            try
                F  = schur(Atmp+complex(eltype(Atmp))(0))
                ps_dict[:schur_mtx] = F.T
                ps_dict[:schur_unitary_mtx] = F.Z
                ps_data.matrix = UpperTriangular(F.T)
                ps_data.unitary_mtx = ps_data.input_unitary_mtx * F.T
                ps_dict[:projection_on] = true
                ps_dict[:ews] = F.values
                ps_dict[:orig_ews] = copy(F.values)
            catch
                ps_data.matrix = Atmp
            end
        end
    else # dense, non-square (but not Hessenberg), and direct
        rfS,rfT = rect_fact(A)
        ps_dict = Dict{Symbol,Any}(:Aisreal => Aisreal,
                                   :isHessenberg => false,
                                   :projection_on => false,
                                   :proj_lev => Inf,
                                   :matrix2 => rfT,
                                   :ews => eigA)
        ps_data = PSAStruct(rfS,input_unitary_mtx,A,input_unitary_mtx,
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
    ps_dict[:threaded] = get(opts,:threaded,false)

    # DEVNOTE: if direct, upstream constructs axes and calls origplot/redraw
    return ps_data
end

"""
    new_matrix(A, opts::Dict{Symbol,Any}=()) -> ps_data

process a linear operator object into the auxiliary data structure used by
Pseudospectra.

There must be methods with `A` for `eltype`, `size`, and `mul!`.
It is up to the user to make sure that `mul!` is consistent with any
options passed to the iterative solver (see documentation for [`xeigs`](@ref)).
"""
function new_matrix(A, opts::Dict{Symbol,Any}=Dict{Symbol,Any}())
    # CHECKME: can A be anything other than a LinearMap here?
    m,n=size(A)
    (m == n) || throw(ArgumentError(
        "Only square linear operators are supported."))
    Aisreal = get(opts,:real_matrix, !(eltype(A) <: Complex))

    verbosity = get(opts,:verbosity,1)
    direct = false
    eigA = get(opts,:eigA,Vector{complex(eltype(A))}())

    input_unitary_mtx = get(opts,:unitary_mtx,I)
    npts = get(opts,:npts,setgridsize(n,24,80,!Aisreal))
    proj_lev = get(opts,:proj_lev,Inf)

    ps_dict = Dict{Symbol,Any}(:Aisreal => Aisreal,
                                   :isHessenberg => false,
                                   :projection_on => false,
                                   :proj_lev => Inf,
                                   :ews => eigA)
    ps_data = PSAStruct(A, input_unitary_mtx, A, input_unitary_mtx,
                        ps_dict)
    if !haskey(opts,:arpack_opts)
        @info("setting default ARPACK options")
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
    ps_dict[:threaded] = false

    return ps_data
end

# for verifying that tests cover intended cases
const logging_algo = Ref{Bool}(false)

"""
    driver!(ps_data, opts, gs=defaultgs(); revise_method=false)

Compute pseudospectra and plot a spectral portrait.

If using an iterative method to get some eigenvalues and a projection, invokes
that first.

# Arguments
- `ps_data::PSAStruct`: ingested matrix, as processed by `new_matrix`
- `gs::GUIState`: object handling graphical output
- `opts::Dict{Symbol,Any}`:
  - `:ax`, axis limits (overrides value stored in `ps_data`).
  - other options passed to `redrawcontour`, `arnoldiplotter!`

Note that many options stored in `ps_data` by `new_matrix()` influence the processing.

When revising a spectral portrait (`revise_method==true`), the following
entries in `opts` also apply:
 - `:arpack_opts::ArpackOptions`,
 - `:direct::Bool`.
"""
function driver!(ps_data::PSAStruct,
                 optsin::Dict{Symbol,Any}=Dict{Symbol,Any}(),
                 gs::GUIState=defaultgs();
                 myprintln=println, logger=:default, revise_method=false)
    # DEVNOTE: mostly corresponds to switch_redraw.m in EigTool
    opts = fillopts(gs,optsin)
    ps_dict = ps_data.ps_dict
    verbosity = get(ps_dict,:verbosity,1)

    # For changing from direct to iterative, or vice versa,
    if revise_method & haskey(opts,:direct)
        set_method!(ps_data, opts[:direct])
    end

    if revise_method & haskey(opts,:arpack_opts) & !ps_dict[:direct]
        ao = opts[:arpack_opts]
        if !isa(ao,ArpackOptions)
            @mywarn(logger,"invalid :arpack_opts option")
            return nothing
        end
        if haskey(ps_dict,:arpack_opts)
            pvalid = (ao == ps_dict[:arpack_opts])
        else
            pvalid = false
        end
        ps_dict[:proj_valid] = pvalid
        ps_dict[:arpack_opts] = ao
    end

    # if caller specifies ax, use it or bust.
    if haskey(opts,:ax)
        if isvalidax(opts[:ax])
            new_ax = opts[:ax]
        else
            @mywarn(logger,"opts[:ax] is not a valid bounding box")
            return nothing
        end
    else
        new_ax = zeros(0)
    end

    if ps_dict[:direct] || get(ps_dict,:proj_valid,false)
        # for iterative methods, we get here on reentrance with the projection
        n,m = size(ps_data.matrix)
        A = ps_data.matrix
        B = get(ps_dict,:matrix2,I)
        eigA = ps_dict[:ews]
        zoom = ps_data.zoom_list[ps_data.zoom_pos]

        if isempty(new_ax)
            if !isempty(eigA)
                # This sets the default domain for the typical case
                isempty(zoom.ax) && (zoom.ax = vec2ax(eigA))
                if !isheadless(gs)
                    # show eigenvalues while waiting
                    ewsplotter(gs, eigA, zoom)
                end
            else
                if isempty(zoom.ax)
                    @mywarn(logger,"bounding box must be specified")
                    return nothing
                end
            end
        else
            zoom.ax = new_ax
        end

        if haskey(opts, :npts)
            new_npts = opts[:npts]
            if isa(new_npts, Integer) && (new_npts > 2) && (new_npts < 2049)
                zoom.npts = new_npts
            else
                @mywarn(logger,"opts[:npts] is not a valid number of points")
                return nothing
            end
        end

        psa_opts = Dict{Symbol,Any}(:levels=>expandlevels(zoom.levels),
                                    :recompute_levels=>zoom.autolev,
                                    :proj_lev=>zoom.proj_lev,
                                    :scale_equal=>zoom.scale_equal,
                                    :real_matrix=>ps_dict[:Aisreal],
                                    :threaded=>ps_dict[:threaded],
                                    :verbosity=>verbosity)
        ss = size(A)
        Z,x,y,t_levels,err,Tproj,eigAproj,algo = psa_compute(A,zoom.npts,
                                                             zoom.ax,
                                                             eigA,psa_opts,
                                                             B,
                                                             logger=logger)

        # FIXME: handle projection properly
        ps_dict[:proj_ews] = eigAproj
        ps_dict[:algo] = algo

        if err != 0
            @mywarn(logger,"PSA computation failed")
            # FIXME: reset GUI if any
            return nothing
        end

        (logging_algo[] | (verbosity > 1)) && println("algorithm: ",algo)
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
        # Iterative method (uses ARPACK):
        # This performs implicitly restart Arnoldi steps to
        # project on a Krylov subspace, yielding a Hessenberg matrix
        # with approximately the same spectral properties (locally) as A.

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
                @warn("xeigs throws:")
                display(JE)
                println()
                stuff = (:failure,nothing)
                put!(chnl,stuff)
                return nothing
            end
        end

        local ews,H,V
        local nconv,niter,nmult
        old_ax = zeros(0)
        (verbosity > 1) && println("calling xeigs w/ $ao")
        chnl = Channel(xeigsproducer)
        xeigs_result = take!(chnl)
        ap_state = nothing
        while xeigs_result[1] ∉ [:finale,:failure]
            the_key,dispvec,the_str,the_ews,the_shifts = xeigs_result
            if !isheadless(gs)
                ap_state = arnoldiplotter!(gs,old_ax,opts,dispvec,
                                           the_str,the_ews, the_shifts,
                                           ap_state)
            end # if gs
            xeigs_result = take!(chnl)
        end
        if xeigs_result[1] == :failure
            @mywarn(logger,"xeigs failed")
            return nothing
        end
        the_key,ews,v,nconv,niter,nmult,resid,H,V = xeigs_result
        if verbosity > 0
            println("xeigs: $nconv of $(ao.nev) converged in $niter iters "
                    * "($nmult MxV)")
        end
        if verbosity > 1
            println("xeigs ews:")
            display(ews); println()
        end
        ews = filter(x->!isnan(x), ews)

        # We basically replace A with H, saving some projection information,
        # and proceed with the dense matrix algorithms.

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
        resize!(ps_data.zoom_list,1)
        # CHECKME: do we need remove() here?
        ps_dict[:mode_markers] = []
        zoom = ps_data.zoom_list[1]

        if isempty(new_ax)
            origax = ps_dict[:init_opts].ax # init_opts is a Portrait!
            if !isempty(origax) && isvalidax(origax) # && !ps_dict[:init_direct]
                copy!(zoom.ax,origax)
            else
                # CHECKME: maybe use init_ews if available?
                zoom.ax = vec2ax(ews)
          # elseif gs.mainph != nothing
          #    println("using eigvals for axis limits")
          #    copy!(zoom.ax,getxylims(gs.mainph))
            end
        else
            zoom.ax = new_ax
        end
        zoom.autolev = ps_dict[:init_opts].autolev
        zoom.levels = deepcopy(ps_dict[:init_opts].levels)

        delete!(ps_dict,:proj_axes)
        delete!(ps_dict,:comp_proj_lev)
        origplot!(ps_data,opts,gs) # WARNING: reenters driver!()
        # caller must reset GUI if appropriate
    end
    nothing
end

function iscomputed(ps_data::PSAStruct, idx=ps_data.zoom_pos)
    ps_data.zoom_list[idx].computed
end

iscomputed(zoom::Portrait) = zoom.computed

"""
Make sure zoom list is ok, then redraw (unless `ax_only`).

Note: truncates zoom list, so use for a new problem or for a reset.
"""
function origplot!(ps_data::PSAStruct, opts, gs; ax_only = false)
    ps_data.zoom_pos = 1
    ps_dict = ps_data.ps_dict
    resize!(ps_data.zoom_list,1)
    zoom = ps_data.zoom_list[1]
    if isempty(zoom.ax)
        if isempty(get(ps_dict,:ews,[]))
            @warn("origplot called w/o preset axes or eigenvalues")
        else
            zoom.ax = vec2ax(ps_dict[:ews])
        end
    end
    ax_only || driver!(ps_data,opts,gs)
    nothing
end

"""
possibly change from direct to iterative method or vice versa
"""
function set_method!(ps_data::PSAStruct, todirect::Bool)
    # this is the part of switch_method which pertains to ps_data
    ps_dict = ps_data.ps_dict
    (todirect == ps_dict[:direct]) && return
    m,n = size(ps_data.input_matrix)
    if todirect
        if haskey(ps_dict,:schur_matrix)
            ps_data.matrix = ps_dict[:schur_matrix]
            ps_data.ews = ps_dict[:orig_ews]
            ps_dict[:ew_estimates] = false
        elseif m==n && issparse(ps_data.input_matrix)
            ps_data.matrix = ps_data.input_matrix
        end
        if haskey(ps_dict,:schur_unitary_mtx)
            ps_data.unitary_mtx  = ps_data.input_unitary_mtx *
                ps_data.schur_unitary_mtx
        else
            ps_data.unitary_mtx = ps_data.input_unitary_mtx
        end
        ps_dict[:proj_valid] = false
        # if reverting to a square matrix, no longer have ARPACK projection
        ss = size(ps_data.matrix)
        (ss[1]==ss[2]) && (ps_dict[:isHessenberg] = false)
    else # switch to iterative
        (m == n) || throw(ArgumentError("Iterative method not implemented "
                                        * "for rectangular matrices"))
        # apparently that's all we need for now
    end
    ps_dict[:direct] = todirect
end

"""
    spectralportrait(A::AbstractMatrix; npts=100) => Plots object

compute pseudospectra of matrix `A` and display as a spectral portrait.

Pseudospectra are computed on a grid of `npts` by `npts` points in
the complex plane, including a neighborhood of the spectrum.
Contour levels are `log10(ϵ)` where `ϵ` is the inverse resolvent norm.
This is a convenience wrapper for simple cases; see the Pseudospectra
package documentation for more elaborate interfaces.
"""
function spectralportrait(A0 :: AbstractMatrix; npts=100)
    pp = getpsplotter()
    local ps_data
    try
        ps_data = new_matrix(A0)
    catch JE
        @warn "The spectralportrait function only works for simple cases."
        rethrow(JE)
    end
    n,m = size(ps_data.matrix)
    A = ps_data.matrix
    ps_dict = ps_data.ps_dict
    B = get(ps_dict,:matrix2,I)
    eigA = ps_dict[:ews]
    if isempty(eigA)
        @error """Unable to proceed without eigenvalues; for non-BLAS types
        import GenericLinearAlgebra and GenericSchur."""
    end
    zoom = ps_data.zoom_list[ps_data.zoom_pos]
    isempty(zoom.ax) && (zoom.ax = vec2ax(eigA))
    psa_opts = _basic_psa_opts(zoom,ps_dict)
    ss = size(A)
    Z,xs,ys,t_levels,err,Tproj,eigAproj,algo = psa_compute(A,npts,
                                                             zoom.ax,
                                                             eigA,psa_opts,
                                                             B)
    return _portrait(pp,xs,ys,Z,eigA)
end

_basic_psa_opts(zoom,ps_dict) = Dict{Symbol,Any}(
    :levels=>expandlevels(zoom.levels),
    :recompute_levels=>zoom.autolev,
    :proj_lev=>zoom.proj_lev,
    :scale_equal=>zoom.scale_equal,
    :real_matrix=>ps_dict[:Aisreal],
    :verbosity=>0)

################################################################
# FIXME: until we think of a better way to handle this:
include("../examples/demo_mtx.jl")

if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base, :get_extension)
function __init__()
    @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("../ext/PseudospectraMakie/PseudospectraMakie.jl")
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("../ext/PseudospectraPlots/PseudospectraPlots.jl")
    @require PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee" include("../ext/PseudospectraPyPlot/PseudospectraPyPlot.jl")
end
end

end # module
