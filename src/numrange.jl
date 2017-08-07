#=
This file is part of Pseudospectra.jl.
See LICENSE file for details.

Julia translation
copyright 2017 Ralph Smith

Portions derived from EigTool
Copyright (c) 2002-2014, The Chancellor, Masters and Scholars
of the University of Oxford, and the EigTool Developers. All rights reserved.
=#

"""
    numrange!(ps_data,npts=20)

compute the numerical range (a.k.a. field of values) of a dense matrix
then store it in the `ps_data`.
"""
function numrange!(ps_data,thmax=20)
    rayleighquotient(B,x) = vecdot(x, B * x) / vecdot(x,x)
    ps_dict = ps_data.ps_dict
    # don't recompute if no change
    if !haskey(ps_dict,:fov) || (length(ps_dict[:fov]) != 2*(thmax+1))
        if !haskey(ps_dict,:schur_mtx)
            throw(ArgumentError("only implemented for Schur-factored matrices"))
        end
        A = ps_dict[:schur_mtx]
        # m,n = size(A)
        T = eltype(A)
        CT = (T <: Real) ? Complex{T} : T
        z = zeros(CT,2*thmax+2)
        for i = 0:thmax
            # upstream allows for interactive cancellation

            # get Hermitian part of rotated A
            th = (i/thmax)*π
            Ath = exp(th*1im)*A
            H = (1/2)*(Ath + Ath')
            d,X = eig(H)

            # RQ's of A correspond to eigenvalues of H w/ extreme real parts
            k = sortperm(real(d))
            z[i+1] = rayleighquotient(A,X[:,k[1]])
            z[1+i+thmax] = rayleighquotient(A,X[:,k[end]])
        end
        z[end] = z[1] # close curve for plotting
        ps_dict[:fov] = z
    else
        z = ps_dict[:fov]
    end

    zoom = ps_data.zoom_list[ps_data.zoom_pos]
    ax = zoom.ax
    if mapreduce(w -> ((real(w) < ax[1]) | (real(w) > ax[2]) |
                       (imag(w) < ax[3]) | (imag(w) > ax[4])),
                 &, true, z)
        warn("The boundary of the numerical range is not visible on the "
             * "current axes; expand axis limits to see it. "
             * "A bounding box is $(extrema(real(z))) $(extrema(imag(z))).")
    end
end
