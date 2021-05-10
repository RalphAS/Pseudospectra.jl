#=
Computational kernels for pseudospectra computations.

This file is part of Pseudospectra.jl.

Julia implementation
Copyright (c) 2017-2021 Ralph A. Smith

Portions derived from EigTool:
 Copyright (c) 2002-2014, The Chancellor, Masters and Scholars
 of the University of Oxford, and the EigTool Developers. All rights reserved.
 EigTool is maintained on GitHub:  https://github.com/eigtool

SPDX-License-Identifier: BSD-3-Clause
License-Filename: LICENSES/BSD-3-Clause_Eigtool
=#

# normally hardwired, but change to get test coverage w/o huge problems
# or to get reference solutions for comparison.
mutable struct ComputeThresholds
    # FIXME: undo mutability when tests etc. are updated
    minlancs4psa::Int # use SVD for n < this
    maxstdqr4hess::Int # use HessQR for n > this (in rectangular case)
    minnev::Int # number of ew's to acquire for projection
    maxit_lancs::Int # bound on Lanczos iterations
end
const _default_thresholds = ComputeThresholds(55,200,20,99)

# FIXME: temporary alias until tests etc. are updated
const psathresholds = _default_thresholds

"""
    psa_compute(T,npts,ax,eigA,opts,S=I) -> (Z,x,y,levels,info,Tproj,eigAproj,algo)

Compute pseudospectra of a (decomposed) matrix.

Uses a modified version of the code in [^Trefethen1999].
If the matrix `T` is upper triangular (e.g. from
a Schur decomposition) the solver is much more efficient than otherwise.

# Arguments
- `T`:      input matrix, usu. from `schur()`
- `npts`:   grid will have `npts × npts` nodes
- `ax`:     axis on which to plot `[min_real, max_real, min_imag, max_imag]`
- `eigA`:   eigenvalues of the matrix, usu. also produced by `schur()`. Pass
  an empty vector if unknown.
- `S`:    2nd matrix, if this is a generalized problem arising from an
  original rectangular matrix.
- `opts`: a `Dict{Symbol,Any}` holding options. Keys used here are as follows:

| Key                 | Type   | Default | Description |
|:-----------|:---------------|:-------------------------------------|:--------|
| `:levels`  | `Vector{Real}` | auto | `log10(ϵ)` for the desired ϵ levels |
| `:recompute_levels` | `Bool` | true | automatically recompute ϵ levels? |
| `:real_matrix`      | `Bool` | `eltype(A)<:Real` | is the original matrix real? (Portrait is symmetric if so.) This is needed because `T` could be complex even if `A` was real.|
| `:proj_lev`         | `Real` |  ∞ | The proportion by which to extend the axes in all directions before projection. If negative, exclude subspace of eigenvalues smaller than inverse fraction. ∞ means no projection.|
| `:scale_equal` | `Bool` | false | force the grid to be isotropic? |
| `:threaded` | `Bool` | false | distribute computation over Julia threads?

# Notes:
- Projection is only done for square, dense matrices.  Projection for sparse
  matrices may be handled (outside this function) by a Krylov method which
  reduces the matrix to a projected Hessenberg form before invoking
  `psa_compute`.
- This function does not compute generalized pseudospectra per se. They may
  be handled by pre- and post-processing.

# Outputs:
- `Z`:         the singular values over the grid
- `x`:         the x coordinates of the grid lines
- `y`:         the y coordinates of the grid lines
- `levels`:   the levels used for the contour plot (if automatically calculated)
- `Tproj`:     the projected matrix (an alias to `T` if no projection was done)
- `eigAproj`:  eigenvalues projected onto
- `algo`: a Symbol indicating which algorithm was used
- `info`:      flag indicating where automatic level creation fails:

| info | Meaning |
|:------|:--------|
| 0 |  No error |
|-1 |  No levels in range specified (either manually, or if matrix is too normal to show levels) |
|-2 |  Matrix is so non-normal that only zero singular values were found |
|-3 |  Computation cancelled |

[^Trefethen1999]: L.N.Trefethen, "Computation of pseudospectra," Acta Numerica 8, 247-295 (1999).
"""
function psa_compute(Targ, npts::Int, ax::Vector, eigA::Vector, opts::Dict, S=I;
                     psatol = 1e-5, thresholds=_default_thresholds,
                     proglog=nothing, logger=:default,
                     ctrlflag=nothing)
    # `proglog` is a placeholder in the hope that something like ProgressMeter can
    #    be made to work in a GUI
    # `ctrlflag` allows user to force early termination.
    m,n = size(Targ)
    eigAproj = copy(eigA) # default
    if isa(S,UniformScaling)
        ms,ns = 1,1
    else
        ms,ns = size(S)
    end
    comp_opts = Dict{Symbol,Any}()
    if !haskey(opts,:recompute_levels)
        comp_opts[:recompute_levels] = false
    end
    if haskey(opts,:levels)
        levels = opts[:levels]
        if length(levels) == 1
            levels = levels * ones(Int,2)
        end
    else
        levels = -8:-1
        comp_opts[:recompute_levels] = true
    end
    all_opts = merge(comp_opts, opts)

    proj_lev = get(all_opts,:proj_lev,Inf)
    re_calc_lev = all_opts[:recompute_levels]
    verbosity = get(all_opts,:verbosity,1)
    threaded = get(all_opts,:threaded,false)

    if get(all_opts,:scale_equal,false)
        y_dist = ax[4]-ax[3]
        x_dist = ax[2]-ax[1]
        if x_dist > y_dist
            x_npts = npts
            y_npts = max(5,ceil(Int,y_dist/x_dist*npts))
        else
            y_npts = npts
            x_npts = max(5,ceil(Int,x_dist/y_dist*npts))
        end
    else
        x_npts = npts
        y_npts = npts
    end
    if get(all_opts,:real_matrix,eltype(Targ)<:Real) && ax[4] > 0 && ax[3] < 0
        y, n_mirror_pts = shift_axes(ax,y_npts)
    else
        n_mirror_pts = 0
        y = collect(range(ax[3], stop=ax[4], length=y_npts))
    end
    x = collect(range(ax[1], stop=ax[2], length=x_npts))
    lx = length(x) # why??
    ly = length(y)
    Z = ones(ly,lx) .+ Inf

    if !issparse(Targ) && n==m
        Tproj, eigAproj = _maybe_project(Targ, proj_lev, ax, eigA, thresholds, verbosity)
        m = size(Tproj,1)
    else
        # sparse or rectangular
        Tproj = Targ
    end # projection branch

    # compute resolvent norms

    local progmeter

    maxit = thresholds.maxit_lancs
    warnflags = falses(2)

    if issparse(Targ)
        algo = :sparse_direct
        Tproj = Targ
        # large value used when subspace eigenproblem doesn't converg
        bigσ = 0.1*floatmax(real(eltype(Targ)))

        if proglog === nothing
            progmeter = Progress(ly,1,"Computing pseudospectra...", 20)
        end
        # reverse order so first row is likely to have a complex gridpt
        # (better timing for LU)
        for j=ly:-1:1
            # check for stop/cancel
            if (ctrlflag !== nothing) && (ctrlflag[] == 1)
                return nothing
                # Eigtool also allows for pause.
            end

            # loop over points in x-direction
            for k=1:lx
                zpt = x[k] + y[j]*im
                t0 = time()
                σ = _psa_lanczos_sparse(Targ, S, zpt, maxit, bigσ)
                Z[j,k] = 1/sqrt(σ)

            end # for k=1:lx
            if proglog === nothing
                update!(progmeter,ly-j+1)
            end
        end
    else # matrix is dense
        step = _get_step_size(m,ly,real(eltype(Tproj)))
        if proglog === nothing
            progmeter = Progress(ly,1,"Computing pseudospectra...", 20)
        end
        for j=ly:-step:1
            # check for stop/cancel
            if (ctrlflag !== nothing) && (ctrlflag[] == 1)
                return nothing
                # Eigtool also allows for pause.
            end

            last_y = max(j-step+1,1)
            q = randn(n) + randn(n)*im
            q = q / norm(q)
            t0 = time()
            Z[j:-1:last_y,:],algo,warnflags = psacore(Tproj,S,q,x,y[j:-1:last_y], m-n+1;
                                                      tol=psatol, threaded=threaded,
                                                      warned=warnflags,
                                                      thresholds=_default_thresholds)

            if proglog === nothing
                update!(progmeter, ly-j+1)
            end
        end # ly loop
    end # if sparse/dense

    # map data (and y) if accounting for symmetry
    if n_mirror_pts < 0
        # bottom half is master
        Z = vcat(Z,reverse(Z[end+n_mirror_pts+1:end,:],dims=1))
        y = vcat(y,-reverse(y[end+n_mirror_pts+1:end]))
    else
        if y[1] != 0
            Z = vcat(reverse(Z[1:n_mirror_pts,:],dims=1),Z)
            y = vcat(-reverse(y[1:n_mirror_pts]),y)
        else
            Z = vcat(reverse(Z[2:n_mirror_pts+1,:],dims=1),Z)
            y = vcat(-reverse(y[2:n_mirror_pts+1]),y)
        end
    end
    ps_tiny = 10*sqrt(floatmin(eltype(Z)))
    (verbosity > 1) && println("range of Z: ",extrema(Z))
    clamp!(Z,ps_tiny,Inf)

    err = 0
    # maybe recalc levels
    if re_calc_lev
        levels,err = recalc_levels(Z,ax)
        if err != 0
            if err == -1
                @mywarn(logger,"Range too small---no contours to plot. Refine grid or zoom out.")
            elseif err == -2
                @mywarn(logger,"Matrix too non-normal---resolvent norm is "
                * "computationally infinite within current axes. Zoom out!")
            end
            return Z,x,y,levels,err,Tproj,eigAproj,algo
        end
    else
        # check that user-supplied levels will plot something
        if ((minimum(levels) > log10(maximum(Z)))
            | (maximum(levels) < log10(minimum(Z))))
            levels, err = recalc_levels(Z,ax)
            @mywarn(logger,"No contours to plot in requested range; 'Smart' levels used.")
            return Z,x,y,levels,err,Tproj,eigAproj,algo
        end
    # check range of Z
        if minimum(levels) < log10(ps_tiny)+1
            @mywarn(logger,"Smallest level allowed by machine precision reached; "
            * "levels may be inaccurate.")
        end
    end
    return Z,x,y,levels,err,Tproj,eigAproj,algo
end

function _get_step_size(n,ly,T)
    nsmall =  precision(T) <= 53 ? 8 : 20
    if n < 100
        step = max(1,floor(Int,ly/nsmall))
    else
        step = min(ly,max(1,floor(Int,4*ly/n)))
    end
    # upstream decreases by factor of 4 if fast implementation is missing
    return step
end

# Trefethen/Wright projection scheme:
# restrict to interesting subspace by ignoring eigenvectors whose
# eigenvalues lie outside rectangle around current axes
function _maybe_project(Targ, factor, ax, eigA, thresholds, verbosity)
    m,n = size(Targ)
    axis_w = ax[2]-ax[1]
    axis_h = ax[4]-ax[3]
    if factor >= 0
        proj_w = axis_w * factor
        proj_h = axis_h * factor
    else
        proj_size = -1 / factor
    end
    np = 0
    ew_range = ax
    # iteratively extend range until 20 (or all) ews are included
    if (m > thresholds.minnev) && !isempty(eigA)
        local selection
        while np < thresholds.minnev
            if factor >= 0
                ew_range = [ew_range[1] - proj_w, ew_range[2] + proj_w,
                            ew_range[3] - proj_h, ew_range[4] + proj_h]
                selection = findall((real(eigA) .> ew_range[1])
                                    .& (real(eigA) .< ew_range[2])
                                    .& (imag(eigA) .> ew_range[3])
                                    .& (imag(eigA) .< ew_range[4]))
            else
                selection = findall(abs.(eigA) .> proj_size)
                proj_size *= (1/2)
            end
            np = length(selection)
            if factor == 0
                # restrict to ews visible in window
                break
            end
        end

    else
        np = m
    end
    # if no need to project (all ews in range)
    if m == np
        m = size(Targ,1)
        # if !opts[:no_waitbar]
        # TODO: post waitbar
        # end
        eigAproj = copy(eigA)
        Tproj = Targ # no mutation, so just dup binding
    else
        if verbosity > 1
            println("projection reduces rank $m -> $np")
        end
        m = np
        n = np
        # restrict eigenvalues and matrix
        eigAproj = eigA[selection]

        # temporarily lose triangular structure
        Tproj = copy(Matrix(Targ))
        # if we have some eigenvalues in our window
        if m>0
            # TODO: post waitbar

            # do the projection
            for i=1:m
                for k=selection[i]-1:-1:i
                    G,r = givens(conj(Tproj[k,k+1]),
                                 conj(Tproj[k,k]-Tproj[k+1,k+1]),
                                 k+1,k)
                    rmul!(Tproj,adjoint(G))
                    lmul!(G,Tproj)
                end
                # TODO: update waitbar
                # TODO: check for pause ll 291ff
                # TODO: check for stop/cancel
            end
            Tproj = UpperTriangular(triu(Tproj[1:m,1:m]))
        end
    end
    return Tproj, eigAproj
end

"""
    psacore(T,S,q,x,y,bw;tol=1e-5,threaded=false) -> Z,algo,warninfo

Compute pseudospectra of a dense triangular matrix

# Arguments
- `T::Matrix{Number}`: long-triangular matrix whose pseudospectra to compute
- `S`: 2nd matrix from generalised pencil `zS-T`. Set to `I` if
           the problem is not generalised
- `q::Vector{Number}`: starting vector for the inverse-Lanczos iteration
           (the same vector is used to start each point in the
           grid defined by `x` and `y`). `q` **must be normalised to
           have unit length.**
- `x::Vector{Real}`: real-part grid to compute the pseudospectra over
- `y::Vector{Real}`: imaginary-part grid to compute the pseudospectra over
- `tol::Real=1e-5`:  tolerance to use to determine when to stop the
           inverse-Lanczos iteration
- `bw::Int`: lower bandwidth of the input matrix (2 for Hessenberg)
- `threaded::Bool`: whether to use multithreading

If `threaded` is `true`, computation of `Z`-values is distributed over
multiple threads. This is worthwhile if using extended precision or a
fine `Z` mesh.  If using BLAS element types, beware of
oversubscription.

# Result
- `Z::Matrix{Real}`: the singular values corresponding to the grid points `x` and `y`.
- `algo::Symbol`: indicates algorithm used
- `warninfo::Vector{Bool}`: records whether warnings were issued
"""
function psacore(T, S, q0, x, y, bw; tol = 1e-5, threaded=false,
                 warned=falses(2), thresholds=_default_thresholds)
    if isreal(T)
        Twork = T .+ complex(eltype(T))(0)
    else
        Twork = copy(T)
    end
    lx = length(x)
    ly = length(y)
    m,n = size(Twork)

    if m<n
        throw(ArgumentError("Matrix size must be m x n with m >= n"))
    end

    generalized =  !isa(S,UniformScaling)
    if generalized
        ms,ns = size(S)
        if (ms != m) || (ns != n)
            throw(ArgumentError("Dimension mismatch for S & T"))
        end
    end

    Z = zeros(ly,lx)
    diaga = diag(Twork)

    # large value used when subspace eigenproblem doesn't converg
    bigσ = 0.1*floatmax(real(eltype(T)))

    # for small matrices just use SVD
    if n < thresholds.minlancs4psa
        Twork = Matrix(Twork)
        if !generalized
            algo = :SVD
            if threaded
                # WARNING: relies on implementation of @threads
                # as in Julia up to (at least) v1.6
                # (stickiness and within-block ordering)
                nt = Threads.nthreads()
                Tw = [similar(Twork) for _ in 1:nt]
                for j=1:ly
                    Threads.@threads for k=1:lx
                        Twt = Tw[Threads.threadid()]
                        zpt = x[k] + y[j]*im
                        copy!(Twt, Twork)
                        Twt[1:m+1:end] .= diaga .- zpt
                        F = svd!(Twt)
                        Z[j,k] = minimum(F.S)
                    end
                end
            else
                for j=1:ly
                    for k=1:lx
                        zpt = x[k] + y[j]*im
                        Twork[1:m+1:end] .= diaga .- zpt
                        F = svd(Twork)
                        Z[j,k] = minimum(F.S)
                    end
                end
            end
        else
            algo = :SVD_gen
            for j=1:ly
                for k=1:lx
                    zpt = x[k] + y[j]*im
                    A = Twork .- zpt*S
                    F = svd!(A)
                    Z[j,k] = minimum(F.S)
                end
            end
        end
    else
        maxit = thresholds.maxit_lancs

        if (m==n) && threaded
            algo = :sq_lanc
            nt = Threads.nthreads()
            Twl = [copy(Twork) for _ in 1:nt]
            Hw = [zeros(real(eltype(Twork)),maxit+1,maxit+1) for _ in 1:nt]
            for j=1:ly
                Threads.@threads for k=1:lx
                    id = Threads.threadid()
                    Twlt = Twl[id]
                    Hwt = Hw[id]
                    zpt = x[k]+y[j]*im
                    Twlt[1:m+1:end] .= diaga .- zpt
                    F1 = UpperTriangular(Twlt)
                    σ, conv1, conv2 = _psa_lanczos!(Hwt, F1, q0, tol, maxit, bigσ)
                    Z[j,k] = 1/sqrt(σ)
                end
            end
        else
          H = zeros(real(eltype(Twork)),maxit+1,maxit+1)
          if m==n
            T1 = copy(Twork)
          end
          for j=1:ly
            for k=1:lx
                zpt = x[k]+y[j]*im
                if m != n
                    if !generalized
                        Twork[1:m+1:end] = diaga .- zpt
                        T1 = copy(Twork)
                    else
                        T1 = Twork - zpt*S
                    end
                    # for large rectangular Hessenberg, use HessQR algorithm
                    if (bw == 2) && (m > thresholds.maxstdqr4hess)
                        algo = :HessQR
                        for jj=1:n-1
                            # DEVNOTE: not using A_mul_B!(G,T1)
                            # because we don't want to mutate top of T1
                            # and views can be expensive
                            G,r = givens(T1[jj,jj],T1[jj+1,jj],1,2)
                            T1[jj:jj+1,jj:end] = G * T1[jj:jj+1,jj:end]
                        end
                    else
                        Qtmp,T1 = qr(T1)
                        algo = generalized ? :rect_qz : :rect_qr
                    end
                    T1 = triu(T1[1:n,1:n])
                else # square
                    algo = :sq_lanc
                    T1[1:m+1:end] .= diaga .- zpt
                end
                if !istriu(T1)
                    F1 = factorize(T1)
                else
                    F1 = UpperTriangular(T1)
                end
                σ, conv1, fail1 = _psa_lanczos!(H, F1, q0, tol, maxit, bigσ)
                if (!conv1) && (!warned[1])
                    @warn "Lanczos convergence failure(s) while computing resolvent norms"
                    warned[1] = true
                end
                if fail1 && (!warned[2])
                    @warn "Eigenvalue convergence failure(s) while computing resolvent norms"
                    warned[2] = true
                end
                Z[j,k] = 1/sqrt(σ)
            end
          end
        end
    end # svd/lanczos branch
    return Z,algo,warned
end

function _psa_lanczos!(H, F1, q0, tol, maxit, bigσ)
    m,n = size(F1)
    # q0 may be too long because of projection
    q = q0[1:n]
    qold = fill!(similar(q),zero(eltype(q)))
    β = 0.0
    σold = 0.0
    local σ
    lancz_converged = false
    H_eigs_failed = false
    for l=1:maxit
        v = (F1 \ (F1' \ q)) - β * qold
        α = real(dot(q,v)) # (q' * v)
        v .= v .- α .* q
        β = norm(v)
        copy!(qold,q)
        q .= v ./ β
        H[l+1,l] = β
        H[l,l+1] = β
        H[l,l] = α
        try
            ewp = eigvals(H[1:l,1:l])
            # eigvals may return complex eltype even if actually real
            σ = maximum(real.(ewp))
            # if !all(isreal.(ewp))
            #     @warn "psa-lancs: eigval anomaly $ewp"
            #     # Should we ask users to report this?
            # end
        catch JE
            # We want a fallback for convergence failure, but throw in
            # other cases.
            # WARNING: this is fragile, depends on library internals
            if !isa(JE, LinearAlgebra.LAPACKException)
                rethrow(JE)
            end
            H_eigs_failed = true
            σ = bigσ
            break
        end
        if (abs(σold / σ - 1) < tol || β == 0)
            lancz_converged = true
            break
        end
        σold = σ
    end
    return σ, lancz_converged, H_eigs_failed
end

function _psa_lanczos_sparse(Targ, S, zpt, maxit, bigσ)
    m,n = size(Targ)
    Tc = complex(eltype(Targ))

    F = lu(Targ - zpt*S)
    σold = 0
    qold = zeros(m)
    β = 0
    H = zeros(Tc,1,1)
    q = normalize!(randn(n) + randn(n)*im)
    w = similar(q)
    v = similar(q)
    local σ
    for l=1:maxit
        ldiv!(w,F,q)
        ldiv!(v,adjoint(F),w)
        v = v - β * qold
        α = real(dot(q,v))
        v = v - α * q
        β = norm(v)
        qold = q
        q = v * (1 / β)
        Hold = H
        H = zeros(Tc,l+1,l+1)
        copyto!(view(H,1:l,1:l),Hold)
        H[l+1,l] = β
        H[l,l+1] = β
        H[l,l] = α
        # calculate eigenvalues of H
        # if error is too big, just set a large value
        try
            ew = eigvals(H[1:l,1:l])
            σ = maximum(real.(ew))
        catch JE
            σ = bigσ
            break
        end
        if (abs(σold / σ - 1)) < 1e-3
            break
        end
        σold = σ
    end
    return σ
end
