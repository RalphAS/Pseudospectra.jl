#=
This file is part of Pseudospectra.jl.

Julia translation
Copyright © 2017 Ralph Smith

Portions from original MATLAB code (via EigTool)
Copyright © 2002-2014 James Burke, Adrian Lewis, Emre Mengi and Michael Overton

SPDX-License-Identifier: BSD-3-Clause
License-Filename: LICENSES/BSD-3-Clause_Eigtool
=#

"""
    psa_abscissa(A,ϵ [,d]) -> α,z

Compute ϵ-pseudospectral abscissa for a dense matrix.

Quadratically convergent two-way method to compute the
ϵ-pseudospectral abscissa `α` of a dense matrix `A`. Also returns a vector
`z` of points where the pseudospectrum reaches the abscissa.
Uses the criss-cross algorithm of Burke et al.

The ϵ-pseudospectral abscissa is
```
   maximum(real(z)) for z s.t. minimum(σ(A-zI)) == ϵ
```

Optional arg:

* `d`: eigenvalues of A, if known in advance

Keyword args:

* `verbosity: 0 for quiet, 1 for noise
"""
function psa_abscissa end

#=
Undocumented keyword args:
* `plotfig`: if a +ve integer, a figure number to use for diagnostic plots (**uses PyPlot**)
* `iterprompt`: whether to pause between steps; for plot visibility

Plotting only works if PyPlot bindings are imported into this module, which we
won't do by default.
=#


function psa_abscissa(A,epsln,eA=[];verbosity=0,plotfig=0,iterprompt=false)
# A = ps_data.input_matrix

# This is the version from EigTool.
#  It is suited to epsilon not too small (otherwise see Overton's web
#  pages, or perhaps Mengi's) and moderate size dense A
#  For large (esp. sparse) A, see Overton's psapsr code.

    if (plotfig > 0)
        if (!isdefined(Main, :PyPlot) || !isdefined(:plot)
            || !(plot === Main.PyPlot.plot))
            @warn("plotting is only implemented for PyPlot")
            plotfig = 0
        end
    end

    n,m = size(A)
    (n==m) || throw(ArgumentError("matrix must be square"))
    (isa(epsln,Real) && (epsln >= 0)) || throw(ArgumentError("ϵ must be >= 0"))

    isempty(eA) && (eA = eigvals(A))

    if plotfig > 0
        figure(plotfig)
        clf()
        plot(real(eA),imag(eA),"k.")
    end

    smalltol = 1e-10 * max(norm(A),epsln)
    purespa = maximum(real.(eA)) # spectral abscissa

    if epsln == 0 # degenerate case
        idx = sortperm(-real.(eA))
        eA = eA[idx]
        f = purespa
        z = [ew for ew in eA if (real(ew) >= (f-smalltol))]
        return f,z
    else
        if (epsln < 1e-9)
            # TODO: divert to specialized Hamiltonian method
            @warn("epsln too small; expect poor accuracy")
        end
        xold = -Inf
        x,ind = findmax(real.(eA)) # initial iterate
        y = imag(eA[ind])
        (verbosity > 0) && @printf("\npsa_abscissa: x=%22.15f   ",x)
        iter = 0
        no_imageig = false
        ybest = NaN
        E = Matrix(epsln * I, size(A)...)
        Areal = (eltype(A) <: Real)
        imagtol = smalltol # used to detect zero real parts

        while !no_imageig && (x > xold)
            iter += 1
            (iter > 20) && error("psa_abscissa: too many steps")
            yold = y
            # given current x, look for all relevant interseections of vertical
            # line w/ pseudospectrum, process and return pair midpoints.
            # note scalar x, scalar ybest, vector y
            # (updates imagtol)
            y,newitol = pspa_2way_imag(A,E,epsln,x,ybest,iter,imagtol,plotfig,verbosity)
            imagtol = newitol
            ptout = "Computing Pseudospectral Abscissa..(iteration $iter)"
            if iterprompt
                print(" [RET]")
                readline()
            end
            no_imageig = isempty(y)
            if !no_imageig
                (verbosity > 0) && print("psa_abscissa: y = $y  ")

                # given resulting y values, look for rightmost intersection of
                # corresponding horizontal lines w/ pseudospectrum
                # note: scalar x, vector y, scalar ybest (even if tied by conj. pair

                xold = x
                ybestt = ybest
                x,ybest = pspa_2way_real(A,E,y,imagtol,plotfig,xold)
                if verbosity > 0
                    println()
                    print("psa_abscissa: x = $x   ")
                end
                if x < xold
                    x = xold # terminates while loop
                    ybest = ybestt
                    if verbosity > 0
                        println()
                        print("psa_abscissa: could not find bigger x")
                    end
                end # if x < xold
            end # if !no_imageig
        end # while

        if isempty(ybest)
            # TODO: fall back to Hamiltonian scheme
            error("Failed in first iteration")
        end

        f = x
        z = x + ybest*im
        if ybest != 0 # FIXME: abs(ybest) < tol
            return f,[z,conj(z)]
        else
            return f,[z]
        end
    end
end

function pspa_2way_real(A,E,y,imagtol,plotfig,xold)

    n = size(A,1)
    Aisreal = (eltype(A) <: Real)
    xnew = zeros(size(y))
    for j in 1:length(y)
        B2 = A - y[j]*im*I
        M2 = vcat(hcat(1im*B2',E),hcat(-E,1im*B2))
        eM2 = eigvals(M2)
        println("test ",minimum(abs.(real.(eM2)))," $imagtol")
        if minimum(abs.(real.(eM2))) <= imagtol # check for real eigenvalue
            xnew[j] = maximum([imag(ew) for ew in eM2 if (abs(real(ew)) < imagtol)]
                              )
        else
            xnew[j] = -Inf
            error("horizontal search failed: no intersection points found for one of the midpoints")
            # TODO: fall back to alternate version
        end
        if plotfig > 0
            figure(plotfig)
            hold(true)
            plot([xold, xnew[j]],y[j]*ones(2),"m-") # horiz line
            plot([xnew[j]],[y[j]],"b+") # right end point
            if Aisreal
                plot([xold,xnew[j]], -y[j]*ones(2),"m-")
                plot([xnew[j]],[-y[j]],"b+")
            end
        end
    end
    xbest, ind = findmax(xnew) # furthest to right
    ybest = y[ind]
    if plotfig > 0
        plot(xbest,ybest,"b*")
    end
    return xbest, ybest
end

"""
 Search for intersections between vertical line with given x component
 and the pseudospectrum.  Start by looking for imaginary eigenvalues of
 Hamiltonian matrix; if there are none, return the empty vector for `ynew`.
 Otherwise, remove any non-extreme imaginary eigenvalues that don't correspond
 to the smallest singular value.  Then sort the eigenvalues into pairs.
 As do so, ensure that `ywant` (a scalar) is in the eigenvalue list, and if
 not, add an extra pair above and below `ywant` (unless the pair containing
 `ywant` is already a very short interval, indicating near convergence to a
 maximizer).  The omission must be caused by rounding, and if we overlook
 this, the process could terminate to a local minimizer of the
 pseudospectrum instead of a global maximizer (consider minus
 the 5,5 Demmel matrix with `ϵ = 0.01`).
 Finally, return the pair midpoints; there must be at least one.

 Updates `imagtol`
"""
function pspa_2way_imag(A,E,epsln,x,ywant,iter,imagtol,plotfig,verbosity)
    svd_check = true
    Areal = (eltype(A) <: Real)
    n = size(A,1)
    B = A - x*I
    # compute eigenvalues of Hamiltonian matrix M = [-B' E; -E B]
    # where E = ϵ I
    M = vcat(hcat(-B',E),hcat(-E, B))
    eM = eigvals(M)
    minreal = minimum(abs.(real.(eM)))
    if (iter == 1) && (minreal > imagtol)
        imagtol += minreal
    end

    newimagtol = imagtol
    ynew = zeros(0)
    if minreal <= imagtol # check if M has an imaginary eigenvalue
        y = sort(imag.([ew for ew in eM if abs(real(ew)) <= imagtol]))
        ny = length(y)
        if svd_check
            indx2 = [1] # check out non-extreme imaginary parts
            for check in 2:(ny-1)
                j = check
                Ashift = A - (x + im*y[j])*I
                u,s,v = svd(Ashift)
                minval,minind = findmin(abs.(s .- epsln))
                if minind == n
                    push!(indx2,j)
                end
            end
            push!(indx2,ny)
            removed = ny - length(indx2)
            if removed > 0
                if verbosity > 0
                    println()
                    print("pspa_2way_imag: singular value test removed $removed eigenvalues")
                end
                y = y[indx2]
            end
        end
        if isodd(length(y))
            # TODO: fallback to better version
            error("odd number of intersection points found by vertical search")
        end
        npairs = floor(Int,length(y) / 2)
        if npairs == 0
            # TODO: fallback to better version
            error("all pairs eliminated by singular value test")
        end

        # organize in pairs and take midpoints
        for j in 1:npairs
            ylow = y[2*j-1]
            yhigh = y[2*j]
            # before taking midpoint, if interval is not very short,
            # check and see if ywant is in this interval, well away from the
            # end points. If so, break pair into two pairs, one above and
            # one below ywant, and take midpoints of both.
            inttol = 0.01 * (yhigh - ylow)
            if (ywant > ylow + inttol) && (ywant < yhigh - inttol)
                push!(ynew,(ylow+ywant)/2)
                push!(ynew,(yhigh+ywant)/2)
            else
                push!(ynew,(ylow+yhigh)/2)
            end
        end
        if plotfig > 0
            figure(plotfig)
            plot(x*ones(2),[maximum(y),minimum(y)],"r-") # vertical line
            plot(x*ones(length(y)),y,"g+") # intersections
            plot(x*ones(length(ynew)),ynew,"bx")
        end

        if Areal
            ynew = [yy for yy in ynew if yy >= 0]
        end
    end
    return ynew, imagtol
end
