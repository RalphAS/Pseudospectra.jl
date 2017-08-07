#=
This file is part of Pseudospectra.jl.
See LICENSE file for details.

Julia translation
copyright 2017 Ralph Smith

Portions from original MATLAB code (via EigTool)
Copyright © 2002-2014 Michael Overton and Emre Mengi
=#

"""
    pspr_2way(A,ϵ [,d]) -> r,z

Compute ϵ-pseudospectral radius for a dense matrix.

Quadratically convergent two-way method to compute the
ϵ-pseudospectral radius `r` of a dense matrix `A`. Also returns a vector
`z` of points where the pseudospectrum intersects the circle of radius `r`.
Uses the "criss-cross" algorithm of Overton and Mengi.

The ϵ-pseudospectral radius is
```
   maximum(abs(z)) for z s.t. minimum(σ(A-zI)) == ϵ
```

Optional arg:

* `d`: eigenvalues of A, if known in advance

Keyword args:

* `verbosity: 0 for quiet, 1 for noise
"""
function pspr_2way end

#=
Undocumented keyword args:
* `plotfig`: if a +ve integer, a figure number to use for diagnostic plots (**uses PyPlot**)
* `iterprompt`: whether to pause between steps; for plot visibility

Plotting only works if PyPlot bindings are imported into this module, which we
won't do by default.
=#

function pspr_2way(A,epsln,eA=zeros(complex(eltype(A)),0);
                   verbosity=0,plotfig=0,iterprompt=false)

# This is the version from EigTool.
#  It is suited to epsilon not too small (otherwise see Overton's web
#  pages, or perhaps Mengi's) and moderate size dense A
#  For large (esp. sparse) A, see Overton's psapsr code

    if (plotfig > 0)
        if (!isdefined(Main, :PyPlot) || !isdefined(:plot)
            || !(plot === Main.PyPlot.plot))
            warn("pspr_2way: plotting is only implemented for PyPlot")
            plotfig = 0
        end
    end

    n,m = size(A)
    (n==m) || throw(ArgumentError("matrix must be square"))
    (isa(epsln,Real) && (epsln >= 0)) || throw(ArgumentError("ϵ must be >= 0"))

    isempty(eA) && (eA = eigfact(A)[:values])

    if plotfig > 0
        figure(plotfig)
        clf()
        plot(real(eA),imag(eA),"k.")
    end

    smalltol = 1e-10 * max(norm(A),epsln)
    purespr = maximum(abs.(eA)) # spectral radius

    if epsln == 0 # degenerate case
        idx = sortperm(-abs.(eA))
        eA = eA[idx]
        f = purespr
        z = [ew for ew in eA if (abs(ew) >= (f-smalltol))]
        return f,z
    else
        if (epsln < 1e-9)
            # TODO: divert to specialized Hamiltonian method
            warn("epsilon too small; expect poor accuracy")
        end
        rold = -1e-6
        r,ind = findmax(abs.(eA)) # initial iterate
        theta = angle.(eA[ind])
        (length(theta) == 1) && (theta = [theta])
        thetaold = theta
        iter = 0
        no_imageig = false
        thetabest = NaN
        mE = epsln * eye(size(A,1))
        Areal = (eltype(A) <: Real)
        realtol = smalltol # used to detect zero real parts
        radtol = smalltol # determines how far ew magnitudes can be from unity

        while !no_imageig && (r > rold)
            if verbosity > 0
                @printf("\npspr_2way: r=%22.15f   ",r)
                @printf("\npspr_2way: theta=")
                print(theta)
            end
            iter += 1
            (iter > 20) && error("pspr_2way: too many steps")
            thetabestt = thetabest
            rold = r
            # given the resulting directions in theta (computed in previous
            # iteration except when iter=1), look for the circle with the
            # greatest radius intersecting pseudo-spectrum boundary.
            # Note: input theta is a vector, but r is scalar
            r,thetabest = pspr_2way_rad(A,mE,theta,realtol,iter,rold,plotfig,verbosity)
            ptout = "Computing Pseudospectral Radius..(iteration $iter)"
            if r > rold
                # given current r, look for all relevant intersections of the
                # circle with radius r and the pseudospectrum, process and
                # return pair midpoints
                thetaold = theta
                theta = pspr_2way_theta(A,mE,epsln,r,thetabest,iter,radtol,
                                        smalltol, plotfig,verbosity)
                if iterprompt
                    print(" [RET]")
                    readline()
                end
                no_thetaeig = isempty(theta)
                no_thetaeig && (rold = r)
            else
                thetabest = thetabestt
                r = rold
            end
        end # while

        if isempty(thetabest)
            # TODO: fall back to Hamiltonian scheme
            error("Failed in first iteration")
        end

        f = r
        z = [r * (cos(thetabest) + 1im*sin(thetabest))]
        if Areal && !isreal(z[1])
            z = [z[1],conj(z[1])]
        end
        return f,z
    end
end

"""
Given a set of angles in theta, in each direction finds the point
on the pseudospectrum boundary with the largest radius, i.e.

  `rnew[j] = max |z| s.t. `σ_min(A - z I) = epsln and angle(z) = theta[j]`

`rbest` is the maximum of the largest radius in any direction

          `rbest = max rnew[j], 1<=j<=length(theta)`

`thetabest` is the angle of a point on the boudary with radius `rbest` in one of
the directions in `theta`, i.e. let

 `rbest = rnew[k]`, for some k, 1<=k<=length(theta),`

then `thetabest = theta[k]`.
"""
function pspr_2way_rad(A,mE,theta,realtol,iter,rold,plotfig,verbosity)

    n = size(A,1)
    Aisreal = (eltype(A) <: Real)
    rnew = zeros(theta)
    for j in 1:length(theta)
        K = vcat(hcat(im*exp(im*theta[j])*A',mE),hcat(-mE,im*exp(-im*theta[j])*A))
        eK = eigfact(K)[:values]

        if minimum(abs.(real.(eK))) <= realtol # check for real eigenvalue
            rnew[j] = maximum([imag(ew) for ew in eK if (abs(real(ew)) < realtol)])
            if plotfig > 0
                figure(plotfig)
                # plot radial line between old and new points
                hold(true)
                plot(rold*cos(theta[j]), rold*sin(theta[j]),"bx")
                plot([rold*cos(theta[j]), rnew[j]*cos(theta[j])],
                     [rold*sin(theta[j]), rnew[j]*sin(theta[j])],"m-")
                plot([rnew*cos(theta[j])], [rnew[j]*sin(theta[j])],"b+")
                if Aisreal
                    plot(rold*cos(theta[j]), -rold*sin(theta[j]),"bx")
                    plot([rold*cos(theta[j]), rnew[j]*cos(theta[j])],
                         [-rold*sin(theta[j]), -rnew[j]*sin(theta[j])],"m-")
                    plot([rnew*cos(theta[j])], [-rnew[j]*sin(theta[j])],"b+")
                end
            end
        else
            # there may be no point on the boundary in this direction
            rnew[j] = -Inf

        end
    end
    if isempty(rnew)  # CHECKME: should this be !isfinite(rnew)
        error("horizontal search failed: no intersection points found for one of the midpoints")
        # TODO: fall back to alternate version
    end
    rbest, ind = findmax(rnew) # outermost
    thetabest = theta[ind]
    if plotfig > 0
        plot(rbest*cos(thetabest),rbest*sin(thetabest),"b*")
    end
    return rbest, thetabest
end

"""
Given a radius `r`, it first computes the intersection points
of eps-pseudospectrum boundary with the circle with radius `r`.
This is achieved by finding the generalized eigenvalues of the
matrix pencil `F - lambda*G` where

       `F=[-eps*I A;r*I 0], G=[0 r*I; A' -eps*I]`

and then performing a singular value test. The singular value test
is neccessary to eliminate the points for which eps is a singular
value but not the smallest one. Finally the midpoints of two
succesive intersection points on the circle with radius r is
calculated, i.e. let `re^(i*theta_i)` and `re^(i*theta_(i+1))` be the
ith and (i+1)th intersection points, then ith midpoint is
re^(i*(theta_i+theta_(i+1))/2). Special attention is paid to
keep the angles in the interval [-pi,pi). Furthermore, as it
was the case for pseudo-absicca code, we specifically considered
the case when the angle from the previous iteration is contained
in one of the intervals. At the exit thetanew contains the
angles of the midpoints.

"""
function pspr_2way_theta(A,mE,epsln,r,thetawant,iter,radtol,smalltol,
                         plotfig,verbosity)
    svd_check = true
    Areal = (eltype(A) <: Real)
    n = size(A,1)

    # compute generalized eigenvalues of matrix pencil F - λ G
    R = r * eye(n)
    O = zeros(n,n)
    F = vcat(hcat(mE,A),hcat(R,O))
    G = vcat(hcat(O,R),hcat(A',mE))
    eM = eigfact(F,G)[:values]

    # extract ews with magnitude (nearly) 1
    eM = [ew for ew in eM if (abs(abs(ew)-1) < radtol)]
    if isempty(eM)
        thetanew = zeros(0)
    else
        idx = sortperm(angle.(eM))
        theta = angle.(eM[idx])
        # Perform singular value test on the points probably on the
        # ϵ-pseudospectrum boundary. The ones actually there have
        # smallest singular value equal to ϵ
        idx2 = Int[]
        for j in 1:length(theta)
            (theta[j] < 0) && (theta[j] += 2π)
            Ashift = A - (r*(cos(theta[j])+im*sin(theta[j]))) * I
            s = svdfact(Ashift)[:S]
            minσ,minind = findmin(abs.(s-epsln))
            (minind == n) && push!(idx2,j)
        end
        nremoved = length(theta) - length(idx2)
        if (nremoved > 0)
            if (verbosity > 0)
                @printf("\npspr_2way_theta: singular value test removed %d ew",nremoved)
            end
            theta = theta[idx2]
        end
        if isempty(theta)
            error("singular value test removed all intersection points: please try a smaller ϵ")
        end

        theta = sort(theta)
        # organize in pairs and take midpoints
        thetanew = zeros(0)
        ind = 0
        # shift thetawant (angle from previous iteration) into canon. interval
        (thetawant < 0) && (thetawant += 2π)

        for j in 1:length(theta)
            thetalow = theta[j]
            if j < length(theta)
                thetahigh = theta[j+1]
            else
                thetahigh = theta[1]+2π
            end
            # before taking the midpoint, if this interval is not very short,
            # check and see if thetawant is in this interval, well away from the
            # end points.  If so, break this pair into two pairs, one above
            # and one below thetawant, and take midpoints of both.
            inttol = 0.01 * (thetahigh - thetalow)

            # needed for last interval
	    if ((thetawant+2π > thetalow + inttol)
                && (thetawant+2*pi < thetahigh - inttol))
	        thetawantt = thetawant + 2π
	    else
	        thetawantt = thetawant
	    end



            if (thetawantt > thetalow + inttol) && (thetawantt < thetahigh - inttol)
                # lower midpoint
	        thetamid = (thetalow + thetawantt)/2

	        # shift thetamid into the interval [-π,π] again
	        (thetamid >= 2π) && (thetamid -= 2π)
	        (thetamid >= π) && (thetamid -= 2π)

	        # remove the midpoint if the minimum singular value is greater
	        # than ϵ, since in this case the midpoint should lie outside the
	        # ϵ-pseudospectrum.
	        if (minimum(svdfact(A-r*exp(im*thetamid)*I)[:S]) <= epsln)
		    push!(thetanew, thetamid)
	        end

                # upper midpoint
                thetamid = (thetawantt + thetahigh)/2

	        # shift thetamid into the interval [-π,π] again
	        (thetamid >= 2π) && (thetamid -= 2π)
	        (thetamid >= π) && (thetamid -= 2π)

	        # remove the midpoint if the minimum singular value is greater
	        # than ϵ, since in this case the midpoint should lie outside the
	        # ϵ-pseudospectrum.
	        if (minimum(svdfact(A-r*exp(im*thetamid)*I)[:S]) <= epsln)
		    push!(thetanew, thetamid)
	        end

            else
                # otherwise, if thetawant is not in the interval
                # take the midpoint of thetalow and thetahigh
                thetamid = (thetalow + thetahigh)/2

	        # shift thetamid into the interval [-π,π] again
	        (thetamid >= 2π) && (thetamid -= 2π)
	        (thetamid >= π) && (thetamid -= 2π)


	        # remove the midpoint if the minimum singular value is greater than
	        # ϵ, since in this case the midpoint should lie outside the
	        # ϵ-pseudospectrum.
	        if (minimum(svdfact(A-r*exp(im*thetamid)*I)[:S]) <= epsln)
		    push!(thetanew, thetamid)
	        end
            end
        end

        if plotfig > 0
            figure(plotfig)
            plot(r*cos(0:0.04:6.4),r*sin(0:0.04:6.4))
            # plot intersection points of circle (radius r) and pseudospectrum
            pointsoncircle = r * (cos.(theta) + im*sin.(theta))
            plot(real(pointsoncircle),imag(pointsoncircle),"ro")
        end

        # DEVNOTE: upstream has this logic but doesn't use ynew
        # Did they mean to replace thetanew?
        #=
        if Areal
            # discard midpoints in lower half plane
            ynew = [yy for yy in thetanew if ((yy >= 0) || (yy == -π))]
        end
        =#
    end
    return thetanew
end
