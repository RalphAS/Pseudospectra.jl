#=
ARPACK wrapper for use in Pseudospectra.jl

This file is part of Pseudospectra.jl, whose LICENSE file applies.

Pseudospectra.jl adaptation
Copyright (c) 2017 Ralph A. Smith

Portions derived from the Base library of Julia:

Copyright (c) 2009-2016: Jeff Bezanson, Stefan Karpinski, Viral B. Shah,
 and other contributors

=#

import Base.LinAlg: BlasInt, ARPACKException, checksquare
import Base.LinAlg.ARPACK: naupd, saupd, eupd_wrapper

"""
    xeigs(A, B, channel=nothing; nev=6, ncv=max(20,2*nev+1), which=:LM,
          tol=0.0, maxiter=300, sigma=nothing, ritzvec=true, v0=zeros((0,)),
          wantH=true, options=nothing, nd=0)

Compute (generalized) eigenpairs of `A x=λ B x` and projections.

Modified version of `eigs()` (q.v.) which optionally
1. provides the projection matrices,
2. provides intermediate states (typically for plotting), and
3. works with a functional form of `A` (not thoroughly tested).

For option (2), pass a `Channel` argument; in this case `xeigs` is
implemented as a producer which fills `channel`.
When finished, it `put`s `(:finale, d,[v,],nconv,niter,nmult,resid[,H,V])` to
`channel`.

For option (3), `A(y,x)` must be a function which applies a linear
operation to `x` and stores the result in `y`, both vectors of length `nd`
(which must be provided as a keyword argument). Only some variants are
implemented for this case.
"""
function xeigs(A, B, chnl = nothing; nd = 0,
               nev::Integer=6, ncv::Integer=max(20,2*nev+1), which=:LM,
               tol=0.0, maxiter::Integer=300, sigma=nothing,
               v0::Vector=zeros((0,)), # zeros(eltype(A),(0,)),
               ritzvec::Bool=true, wantH::Bool=true,
               options = Dict{Symbol,Any}())

    if isa(A,AbstractMatrix)
        n = checksquare(A)

        T = eltype(A)
        iscmplx = T <: Complex
        sym = issymmetric(A) && issymmetric(B) && !iscmplx
    else
        (nd == 0) && throw(ArgumentError("caller must provide nd for function version"))
        n = nd
        T = eltype(v0)
        iscmplx = true
        sym = false
    end
    isgeneral = B !== I
    nevmax=sym ? n-1 : n-2
    if nevmax <= 0
        throw(ArgumentError("Input matrix A is too small. Use eigfact instead."))
    end
    if nev > nevmax
        warn("Adjusting nev from $nev to $nevmax")
        nev = nevmax
    end
    if nev <= 0
        throw(ArgumentError("requested number of eigenvalues (nev) must be ≥ 1, got $nev"))
    end
    ncvmin = nev + (sym ? 1 : 2)
    if ncv < ncvmin
        warn("Adjusting ncv from $ncv to $ncvmin")
        ncv = ncvmin
    end
    ncv = BlasInt(min(ncv, n))
    bmat = isgeneral ? "G" : "I"
    isshift = sigma !== nothing

    if isa(which,AbstractString)
        warn("Use symbols instead of strings for specifying which eigenvalues to compute")
        which=Symbol(which)
    end
    if (which != :LM && which != :SM && which != :LR && which != :SR &&
        which != :LI && which != :SI && which != :BE)
        throw(ArgumentError("which must be :LM, :SM, :LR, :SR, :LI, :SI, or :BE, got $(repr(which))"))
    end
    if which == :BE && !sym
        throw(ArgumentError("which=:BE only possible for real symmetric problem"))
    end
    isshift && which == :SM && warn("use of :SM in shift-and-invert mode is not recommended, use :LM to find eigenvalues closest to sigma")

    if which==:SM && !isshift # transform into shift-and-invert method with sigma = 0
        isshift=true
        sigma=zero(T)
        which=:LM
    end

    if sigma !== nothing && !iscmplx && isa(sigma,Complex)
        throw(ArgumentError("complex shifts for real problems are not yet supported"))
    end
    sigma = isshift ? convert(T,sigma) : zero(T)

    if !isempty(v0)
        if length(v0) != n
            throw(DimensionMismatch())
        end
        if eltype(v0) != T
            throw(ArgumentError("starting vector must have element type $T, got $(eltype(v0))"))
        end
    end

    whichstr = "LM"
    if which == :BE
        whichstr = "BE"
    end
    if which == :LR
        whichstr = (!sym ? "LR" : "LA")
    end
    if which == :SR
        whichstr = (!sym ? "SR" : "SA")
    end
    if which == :LI
        if !sym
            whichstr = "LI"
        else
            throw(ArgumentError("largest imaginary is meaningless for symmetric eigenvalue problems"))
        end
    end
    if which == :SI
        if !sym
            whichstr = "SI"
        else
            throw(ArgumentError("smallest imaginary is meaningless for symmetric eigenvalue problems"))
        end
    end

    # Refer to ex-*.doc files in ARPACK/DOCUMENTS for calling sequence
    if isa(A,Function)
        matvecA! = A
        if isshift
            throw(ArgumentError("functional version only implemented for regular mode"))
        end
    else
        matvecA!(y, x) = A_mul_B!(y, A, x)
    end

    if !isgeneral           # Standard problem
        matvecB = x -> x
        if !isshift         #    Regular mode
            mode       = 1
            solveSI = x->x
        else                #    Shift-invert mode
            mode       = 3
            F = factorize(A - UniformScaling(sigma))
            solveSI = x -> F \ x
        end
    else                    # Generalized eigenproblem
        matvecB = x -> B * x
        if !isshift         #    Regular inverse mode
            mode       = 2
            F = factorize(B)
            solveSI = x -> F \ x
        else                #    Shift-invert mode
            mode       = 3
            F = factorize(A - sigma*B)
            solveSI = x -> F \ x
        end
    end

    # Compute the Ritz values and Ritz vectors
    (resid, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, TOL) =
       xaupd_wrapper(T, matvecA!, matvecB, solveSI, n, sym, iscmplx, bmat, nev, ncv, whichstr, tol, maxiter, mode, v0, options, chnl)

    V = copy(v)

    # Postprocessing to get eigenvalues and eigenvectors
    output = eupd_wrapper(T, n, sym, iscmplx, bmat, nev, whichstr, ritzvec, TOL,
                                 resid, ncv, v, ldv, sigma, iparam, ipntr, workd, workl, lworkl, rwork)

    # Issue 10495, 10701: Check that all eigenvalues are converged
    nev = length(output[1])
    nconv = output[ritzvec ? 3 : 2]
    nev ≤ nconv || warn("not all wanted Ritz pairs converged. Requested: $nev, converged: $nconv")

    if wantH
        H = extract_mtx(matvecA!, V,ipntr,workl,ncv,iscmplx,sym)
        output = (output...,H,V)
    end
    output = (:finale, output...)

    if chnl == nothing
        return (output[2:end]...)
    else
        if VERSION < v"0.6-"
            produce(output)
        else
            put!(chnl,output)
        end
    end
end # xeigs

function extract_mtx(matvecA!, V,ipntr,workl,ncv,iscmplx,sym)
    ih = ipntr[5]

    if !iscmplx
        H = zeros(ncv,ncv)
        if sym
            H = (diagm(workl[ih+1:ih-1+ncv],1)
                 + diagm(workl[ih+ncv:ih-1+2*ncv],0)
                 + diagm(workl[ih+1:ih-1+ncv],-1))
        else
            H[:] = workl[ih:ih-1+ncv^2]
        end
    else
        H = zeros(ncv,ncv) + 0im
        # CHECKME: workl is complex, but does ARPACK set ih for a real array?
        H[:] = workl[ih:ih-1+ncv^2]
    end


    if !iscmplx && sym
        β = workl[ih]
    else
        w = zeros(size(V,1))
        matvecA!(w,V[:,end])
        h = (w' * V)'
        mtxprod = w - V*h
        β = norm(mtxprod)
    end

    # form a rectangular Hessenberg matrix
    H = vcat(H,hcat(zeros(1,ncv-1),β))
end

"""
modified version of aupd_wrapper with optional plotting of intermediate results.
"""
function xaupd_wrapper(T, matvecA!::Function, matvecB::Function, solveSI::Function, n::Integer,
                      sym::Bool, cmplx::Bool, bmat::String,
                      nev::Integer, ncv::Integer, which::String,
                       tol::Real, maxiter::Integer, mode::Integer, v0::Vector,
                       options::Dict,chnl)

    lworkl = cmplx ? ncv * (3*ncv + 5) : (sym ? ncv * (ncv + 8) :  ncv * (3*ncv + 6) )
    TR = cmplx ? T.types[1] : T
    TOL = Array{TR}(1)
    TOL[1] = tol

    v      = Array{T}(n, ncv)
    workd  = Array{T}(3*n)
    workl  = Array{T}(lworkl)
    rwork  = cmplx ? Array{TR}(ncv) : Array{TR}(0)

    if isempty(v0)
        resid  = Array{T}(n)
        info   = zeros(BlasInt, 1)
    else
        resid  = deepcopy(v0)
        info   = ones(BlasInt, 1)
    end
    iparam = zeros(BlasInt, 11)
    ipntr  = zeros(BlasInt, (sym && !cmplx) ? 11 : 14)
    ido    = zeros(BlasInt, 1)

    iparam[1] = BlasInt(1)       # ishifts
    iparam[3] = BlasInt(maxiter) # maxiter
    iparam[7] = BlasInt(mode)    # mode

    zernm1 = 0:(n-1)

    all_shifts = []
    done_one_already = false
    H = []
    V = []
    iterx = 0

    while true
        if cmplx
            naupd(ido, bmat, n, which, nev, TOL, resid, ncv, v, n,
                  iparam, ipntr, workd, workl, lworkl, rwork, info)
        elseif sym
            saupd(ido, bmat, n, which, nev, TOL, resid, ncv, v, n,
                  iparam, ipntr, workd, workl, lworkl, info)
        else
            naupd(ido, bmat, n, which, nev, TOL, resid, ncv, v, n,
                  iparam, ipntr, workd, workl, lworkl, info)
        end
        if info[1] == 1
            warn("incomplete convergence. Try a different"
                 * " starting vector or increase maxiter or ncv.")
            break
#            throw(ARPACKException("incomplete convergence. Try a different"
#                        * " starting vector or increase maxiter or ncv."))
        elseif info[1] != 0
            throw(ARPACKException(info[1]))
        end

        x = view(workd, ipntr[1]+zernm1)
        y = view(workd, ipntr[2]+zernm1)
        if mode == 1  # corresponds to dsdrv1, dndrv1 or zndrv1
            if ido[1] == 1
                matvecA!(y, x)
            elseif ido[1] == 99
                break
            else
                throw(ARPACKException("unexpected behavior"))
            end
        elseif mode == 3 && bmat == "I" # corresponds to dsdrv2, dndrv2 or zndrv2
            if ido[1] == -1 || ido[1] == 1
                y[:] = solveSI(x)
            elseif ido[1] == 99
                break
            else
                throw(ARPACKException("unexpected behavior"))
            end
        elseif mode == 2 # corresponds to dsdrv3, dndrv3 or zndrv3
            if ido[1] == -1 || ido[1] == 1
                matvecA!(y, x)
                if sym
                    x[:] = y    # overwrite as per Remark 5 in dsaupd.f
                end
                y[:] = solveSI(y)
            elseif ido[1] == 2
                y[:] = matvecB(x)
            elseif ido[1] == 99
                break
            else
                throw(ARPACKException("unexpected behavior"))
            end
        elseif mode == 3 && bmat == "G" # corresponds to dsdrv4, dndrv4 or zndrv4
            if ido[1] == -1
                y[:] = solveSI(matvecB(x))
            elseif  ido[1] == 1
                y[:] = solveSI(view(workd,ipntr[3]+zernm1))
            elseif ido[1] == 2
                y[:] = matvecB(x)
            elseif ido[1] == 99
                break
            else
                throw(ARPACKException("unexpected behavior"))
            end
        else
            throw(ArgumentError("ARPACK mode ($mode) not yet supported"))
        end

        # This is just a call counter since
        # standard ARPACK doesn't expose the real iteration nr
        iterx += 1
        # TODO: displayRitzValues

        if done_one_already
            # arpackgui_update
            if !cmplx
                if sym
                    dispvec = workl[ipntr[6]+(0:ncv-1)]
                    if which == :BE
                        # roughly nev large and small ews
                        the_ews = dispvec[end-2*nev+1:end]
                        the_shifts = dispvec[1:end-2*nev]
                    else
                        # nev eigenvalues
                        the_ews = dispvec[end-nev+1:end]
                        the_shifts = dispvec[1:end-nev]
                    end
                else # non-sym real
                # println("ew/shifts:",workl[ipntr[6]+(0:2*ncv-1)])
                    dispvec = (workl[ipntr[6]+(0:ncv-1)]
                               + workl[ipntr[7]+(0:ncv-1)]*im)
#                               + workl[ipntr[6]+(ncv:2*ncv-1)]*im)
                end
                # nev+1 ews (keep conjugate pairs together)
                the_ews = dispvec[end-nev:end]
                the_shifts = dispvec[1:end-nev-1]
            else # complex
                dispvec = workl[ipntr[6]+(0:ncv-1)]
                the_ews = dispvec[end-nev+1:end]
                the_shifts = dispvec[1:end-nev]
            end
            if isempty(all_shifts)
                all_shifts = copy(the_shifts)
            else
                append!(all_shifts,the_shifts)
            end
            # we want to use the previous H for PSA/RVec computation
            # since the current one is computed after the implicit restart
            # while the R vals and shifts apply before it.
            prevH = H
            prevV = V
            isempty(V) && (V = copy(v))
            H = extract_mtx(matvecA!, V,ipntr,workl,ncv,cmplx,sym)
            V = copy(v)

            # DEVNOTE: upstream has logic for analyzing a Ritz value here
            # if user invoked "pause"

            if ((chnl != nothing) && get(options,:ARPACK_plots,true)
                && (get(options,:ProgRVals,true) || get(options,:ProgPSA,false))
                && (mod(iterx, get(options,:ARPACK_plotinterval,ncv-nev)) == 0))
                # we only want to update about once per restart, but ARPACK
                # doesn't report restarts, so default to approximate
                # interval (ncv-nev)

                # DEVNOTE: upstream plots Pseudospectra
                # here if options[:ProgPSA] is set

                the_str = @sprintf("dim = %d, iter %d",n,iterx/(ncv-nev))
                stuff = (:media,dispvec,the_str,the_ews,
                         get(options,:ProgAllShifts,false) ? all_shifts
                                                           : the_shifts)
                if VERSION < v"0.6-"
                    produce(stuff)
                else
                    put!(chnl,stuff)
                end
            end
        end #arpackgui_update
        done_one_already = true

    end # iteration loop

    return (resid, v, n, iparam, ipntr, workd, workl, lworkl, rwork, TOL)
end
