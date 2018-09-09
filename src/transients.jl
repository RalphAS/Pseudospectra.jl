#=
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

"""
    transient_bestlb(ps_data, ftype, it_range) => pts,bnds,sel_pt

compute lower bound for transient growth of a linear system

`ftype` may be `:exp` or `:powers` to indicate the sort of transient.

This is based on Trefethen & Embree, Spectra and Pseudospectra, chs. 15, 16.
"""
function transient_bestlb(ps_data::PSAStruct,ftype::Symbol,it_range;
                          dk=0.25, method=:new)
    A = ps_data.matrix
    zoom = ps_data.zoom_list[ps_data.zoom_pos]
    x,y,Z = zoom.x, zoom.y, zoom.Z
    badresult = [NaN],[NaN],[NaN]
    if ftype == :exp
        R = 1.0 ./ Z # resolvent norm
        ωA = 0.5 * maximum(eigvals(A+A')) # numerical abscissa, Eq. 14.2seq.
        if method == :old
            xcrit = 0
        else
            xcrit = ωA
        end
        xtmp = [ifelse(xi<=xcrit,NaN,xi) for xi in x]
        if all(isnan.(xtmp))
            @warn("No points in portrait reach threshold [Re(z)>" *
                 "$xcrit]: expand grid and recompute")
            return badresult
        end
        X = xtmp .* ones(1,length(y))
        t_pts = collect(it_range)
        nt = length(t_pts)
        bnd = zeros(nt)
        sel_pt = Matrix{eltype(x)}(undef,nt,2)
        for it in eachindex(t_pts)
            t = t_pts[it]
            eat = exp.(t.*X)
            if method == :old
                bnds = eat ./ (1 .+ (eat .- 1) ./ (X .* R)) # Eq. 15.11
            elseif method == :new
                # Eq. 15.20, using K as defined at top of Thm. 15.5
                bnds = eat .- (eat .- exp(ωA*t)) ./ ((X .- ωA) .* R)
            else
                # combine the two
                bnds1 = eat ./ (1 .+ (eat .- 1) ./ (X .* R))
                bnds2 = eat .- (eat .- exp(ωA*t)) ./ ((X .- ωA) .* R)
                bnds = [max(x1,x2) for (x1,x2) in zip(bnds1,bnds2)]
            end
            if all(isnan.(bnds))
                @warn("no meaningful results obtained")
                return badresult
            else
                bnds = map(x->isnan(x) ? 0 : x,bnds)
                m,idx = findmax(vec(bnds))
            end
            j,i = map(x->x+1,divrem(idx-1,size(X,1)))
            bnd[it] = m
            sel_pt[it,:] = [x[i],y[j]]
        end
        data_pts = t_pts
    elseif ftype == :powers
        R = 1.0 ./ Z
        X = x .* ones(1,length(y))
        Y = ones(length(x),1) .* y'
        pts = abs.(X .+ 1im * Y)
        anrm = norm(A)
        if method == :old
            rcrit = 1
        else
            rcrit = anrm
        end
        pts = map(x->ifelse(x<=rcrit,NaN,x), pts)
        if all(isnan.(pts))
            @warn("No points achieve threshold [|z| > $rcrit]; "
                 * "expand grid and recompute")
            return badresult
        end
        k_pts = collect(it_range[1]:dk:it_range[end]) # collect(it_range)
        nt = length(k_pts)
        bnd = zeros(nt)
        sel_pt = Matrix{eltype(x)}(nt,2)
        nA = norm(A)
        for it in eachindex(k_pts)
            k = k_pts[it]
            pk = pts .^k
            if method == :old
                bnds = pk ./ (1+(pk-1)./((pts-1) .* (pts .* R -1))) # Eq. 16.15
            elseif method == :new
                bnds = pk .- (pk .- nA^k) ./ ((pts -nA).*R) # Eq. 16.24
            else
                bnds1 = pk ./ (1+(pk-1)./((pts-1) .* (pts .* R -1))) # Eq. 16.15
                bnds2 = pk .- (pk .- nA^k) ./ ((pts -nA).*R) # Eq. 16.24
                bnds = [max(x1,x2) for (x1,x2) in zip(bnds1,bnds2)]
            end
            m,idx = findmax(vec(bnds))
            j,i = map(x->x+1,divrem(idx-1,size(X,1)))
            bnd[it] = m
            sel_pt[it,:] =  [x[i],y[j]]
        end
        data_pts = k_pts
    else
        throw(ArgumentError("incomprehensible ftype=$ftype"))
    end
    data_pts,bnd,sel_pt
end
