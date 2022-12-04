#=
This file is part of Pseudospectra.jl., whose LICENSE file applies.

Julia translation
copyright (c) 2017 Ralph Smith

Portions derived from EigTool
Copyright (c) 2002-2014, The Chancellor, Masters and Scholars
of the University of Oxford, and the EigTool Developers. All rights reserved.

SPDX-License-Identifier: BSD-3-Clause
License-Filename: LICENSES/BSD-3-Clause_Eigtool
=#

"""
    modeplot(ps_data, pkey [,z])

Extract and plot an eigenmode or pseudo-eigenmode for the matrix
encapsulated in the Pseudospectra object `ps_data`.
Use the value `z` if provided or prompt for one.
If `pkey` is 1, find the pseudoeigenmode for `z`;
otherwise find the eigenmode for the eigenvalue closest to `z`.
"""
function modeplot(ps_data::PSAStruct, pkey, z=NaN; gs=defaultgs(),
                  zgetter=replzdlg, zmarker=addmark, verbosity = 0)
    # based on switch_pseudomode()
    # FIXME: a lot of logic is missing here
    ps_dict = ps_data.ps_dict
    pseudo = (pkey == 1) # in lieu of Bool arg FSO GUI
    # get a point
    if isnan(z)
        z = zgetter(gs, what=(pseudo ? "a z-value" : "approx. eigenvalue"))
    end
    if isnan(z)
        warn("modeplot: invalid z")
        return
    end
    z = complex(z)
    if !pseudo
        eigA = get(ps_dict,:proj_ews,ps_dict[:ews])
        isempty(eigA) && throw(ArgumentError("ps_data holds no eigenvalues"))
        dists = [abs(eig - z) for eig in eigA]
        idx = argmin(dists)
        z = eigA[idx]
    end
    if zmarker !== nothing
        zmarker(gs,z,(pseudo ? :pseudo : :eigen))
    end
    # FIXME: check for direct method, etc.
    if !isempty(get(ps_dict,:schur_mtx,zeros(0)))
        A = ps_dict[:schur_mtx]
        U = ps_dict[:schur_unitary_mtx]
    elseif !isempty(get(ps_dict,:proj_matrix,zeros(0))) # && etc.
            A = ps_dict[:proj_matrix]
            U = ps_dict[:proj_unitary_mtx]
    # elseif sparse_direct (no Schur)
    else
        A = ps_data.input_matrix
        U = ps_data.input_unitary_mtx
    end
    # in Arnoldi case extract upper square mtx
    mm,nn = size(A)
    if mm == nn
        plotmode(gs,z,A,U,pseudo,false,verbosity) # avoid pointless copy
    else
        plotmode(gs,z,A[1:nn,1:nn],U,pseudo,false,verbosity)
    end
    # TODO: store marker info in gs
    nothing
end

"""
    oneeigcond(T,ew,verbosity;dlg=replqdlg) -> cond,evr,evl

Obtain condition and eigenvector(s) for eigenvalue `ew` associated with
Schur factor `T`.
"""
function oneeigcond(T0,ew0,verbosity; dlg=replqdlg, max_its=3)
    quest_str = ("The eigenvector could not be determined by inverse iteration."
                 * "\nDo you want to do a (possibly slow) dense SVD instead?")
    n,m = size(T0)
    @assert n==m
    tol = 1e-14
    # perturb eigenvalue to avoid exact singularity
    ew = ifelse(n<200,ew0,ew0*(1+2*eps(abs(ew0))))
    T = copy(T0)
    T = T - ew * I
    # for small matrices, use SVD
    if n < 200
        U,S,V = svd(Matrix(T))
        evr = V[:,end]
        evl = conj(U[:,end])
    else
        # inverse iteration
        if issparse(T)
            F = lu(T)
        end
        wb_size = 200
        v0 = normalize!(randn(n) .+ 0im)
        infs_found = false # indicator for inv-iter success
        evr = copy(v0)
        local resid
        for i=1:max_its
            oldev = copy(evr)
            try
                if issparse(T)
                    evr = F \ oldev
                else
                    evr = T \ oldev
                end
            catch JE
                # FIXME: check exception type: LAPACKException or ??
                infs_found = true
            end
            any(isinf.(evr)) && (infs_found = true)
            if infs_found
                break
            end
            λ = (oldev' * evr)[1]
            resid = norm(evr - λ * oldev) / abs(λ)
            normalize!(evr)
            (verbosity > 2) && println("evr resid: $resid")
            if resid < tol
                break
            end
        end
        if !(infs_found || (resid > tol))
            # no problem w/ right ev, so do left one
            evl = v0
            for i=1:max_its
                oldev = copy(evl)
                try
                    if issparse(T)
                        ldiv!(evl,transpose(F),oldev)
                    else
                        evl = transpose(T) \ oldev # At_ldiv_B!(evl,T,oldev)
                    end
                catch JE
                    infs_found = true
                end
                any(isinf.(evl)) && (infs_found = true)
                if infs_found
                    break
                end
                λ = (oldev' * evl)[1]
                resid = norm(evl - λ * oldev) / abs(λ)
                normalize!(evl)
                (verbosity > 2) && println("evl resid: $resid")
                if resid < tol
                    break
                end
            end
        end
        if infs_found || (resid > tol)
            (verbosity > 1) && println("Infs found: $infs_found resid: $resid")
            if dlg("Compute using full SVD?",quest_str) == 1 # 1:yes
                U,S,V = svd(Matrix(T))
                evr = V[:,end]
                evl = conj(U[:,end])
            else
                evr = NaN
                evl = NaN
            end
        end
    end
    z0 = transpose(evl)*evr
    the_cond = abs(1/z0[1])
    return the_cond, evr, evl
end

"""
     psmode_inv_lanczos(A,q,z,tol,num_its) -> Z,q

use inverse-Lanczos iteration to compute a pseudomode

 A         the matrix to use to find min σ(zI-A)
 q         starting vector for the iteration
 z         point in the complex plane to use
 tol       tolerance to use
 num_its   Maximum number of iterations

 Z         the smallest singular values
 q         the smallest singular vector
"""
function psmode_inv_lanczos(A,q0::Vector,z::Complex,tol,num_its;
                            fallbacksvd=false)
    m,n = size(A)
    (m == n) || throw(ArgumentError("Matrix must be square"))
    diaga = diag(A)
    cdiaga = conj(diaga)
    q = q0 .+ 0im
    if n < 200
        # just use plain SVD
        A1 = copy(A) .+ 0im
        A1[1:m+1:end] = diaga .- z
        U,S,V = svd(Matrix(A1))
        Z = minimum(S)
        @assert S[end] == Z
        q = V[:,end]
    else
        H = zeros(num_its+1,num_its+1)
        if issparse(A)
            F = lu(A-z*I)
            # use UMFPACK factor method to avoid this:
            # Fp = lu(A'-conj(z)*I)
        else
            T1 = A - z * I
            T2 = T1'
        end
        qold = zeros(n)
        Q = copy(q0)
        β = 0.0
        sigold = 0.0
        local σ,α
        l = 0
        for it=1:num_its
            l += 1
            if issparse(A)
                # v = F \ (Fp \ q) - β * qold
                w = zeros(q+0im)
                v = similar(w)
                ldiv!(w,adjoint(F),q)
                ldiv!(v,F,q)
                v -= β * qold
            else
                v = T1 \ (T2 \ q) - β * qold
            end
            α = real((q' * v)[1])
            v = v - α * q
            β = norm(v)
            qold = q
            q = v / β
            Q = hcat(Q,q)
            H[l+1,l] = β
            H[l,l+1] = β
            H[l,l] = α
            # calculate eigenvalues of H, but if error is too big,
            # set σ to a huge value
            try
                E = eigen(H[1:l,1:l])
                σ = maximum(E.values)
            catch
                σ = 1e308
            end
            resid = abs(sigold/σ - 1)
            if (resid < tol) || (β == 0)
                break
            end
            sigold = σ
        end
        H = H[1:l,1:l]
        Z = 1 / sqrt(σ)
        # get eigenvector from inverse Lanczos basis
        # if matrix is too non-normal, revert
        try
            E = eigen(H)
            q = Q[:,1:end-1] * E.vectors[:,end]
        catch
            l = num_its
        end
        if l >= num_its
            if fallbacksvd
                A1 = A - z * I
                U,S,V = svd(Matrix(A1))
                Z = minimum(S)
                q = V[:,end]
            else
                Z = NaN
                q = NaN
            end
        end
    end # if n < 200
    return Z,q
end
