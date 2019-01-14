#=
Toeplitz matrix examples from Trefethen & Embree, Spectra & Pseudospectra
mainly ch. 7

Instructions
============

1. Establish a plotter for Pseudospectra:
```
using Pseudospectra
setpsplotter(:PyPlot) -or- setpsplotter(:Plots)
```
If using Plots, you should probably set a backend first.

2. `include()` this file.
3. Invoke `M.runme()` with arguments of interest.
=#
module M

using Pseudospectra, ToeplitzMatrices, Polynomials, LinearAlgebra

p = getpsplotter()
if p == :undef
    error("plotter must be set first")
elseif p == :PyPlot
    using PyPlot
else
    using Plots
end

"""
Draws a spectral portrait of a Toeplitz matrix, also showing
the locations of eigenvalues of perturbations of said matrix.
"""
function runme(casename="limaçon",n=64,niter=100)

    plotter = getpsplotter()

    gs = setgs()

    ax = zeros(0) # override as needed

    # conventions: x^0 terms must be equal
    # make life easier for those of us w/o French keyboards
    if casename[1:3] == "lim" # limaçon
        ax = [-2,2,-2,2]
        vc = [0.0,1.0,1.0]
        vr = [0.0,0.0]
    elseif casename[1:3] == "tri" # triangle
        vc = [0.0,1.0]
        vr = [0.0,0.0,1/4]
    elseif casename[1:3] == "wha" # whale
        vr = [0.0,1.0,1im,-3-2im,-1]
        vc = [0.0,10,3+1im,4,1im]
    elseif casename[1:3] == "but" # butterfly
        vc = [0.0,-1im,1.0]
        vr = [0.0,1im,-1]
    elseif casename[1:3] == "grc" # Grcar
        vc = ones(4)
        vr = [1.0,-1.0]
    elseif casename[1:3] == "ch3" # chapter 3 example
        ax = [-1.5,1.5,-1,1]
        vc = [0.0,1/4]
        vr = [0.0,1.0]
    elseif casename[1:3] == "ch7" # chapter 7 example
        vc = [0.0,0,-4,-2im]
        vr = [0.0,2im,-1,2]
    else
        throw(ArgumentError("unrecognized casename"))
    end

    vrp = Poly(vr)
    vcp = Poly(vcat([0.0],vc[2:end]))

    symbolfunc(z) = polyval(vrp,z) + polyval(vcp,1/z)

    At = Toeplitz(vcat(vc,zeros(n-length(vc))),vcat(vr,zeros(n-length(vr))))
    # some methods are missing for Toeplitz, so ...
    A = Matrix(At)

    opts=Dict{Symbol,Any}()
    isempty(ax) || (opts[:ax] = ax)
    ps_data = new_matrix(A,opts)
    driver!(ps_data,opts,gs)

    ϵ = 0.001
    for it=1:niter
        pert = randn(size(A)) + im*randn(size(A))
        pn = norm(pert)
        pert .*= ϵ/pn
        A1 = A + pert
        ews = eigvals(A1)
        if plotter == :PyPlot
            plot(real(ews),imag(ews),"r.")
        else
            scatter!(real(ews),imag(ews),markersize=2,c=:red,leg=false,
                     markerstrokealpha=0.0)
        end
    end
    θv = range(0, stop=2π, length=201)

    symcurve = [symbolfunc(cis(θ)) for θ in θv]
    if plotter == :PyPlot
        plot(real(symcurve),imag(symcurve),"k")
        p = nothing
    else
        p = plot!(real(symcurve),imag(symcurve))
    end

    #=
    # Sometimes perturbed ews are inside this curve.
    # Theory seems confused.

    efac = ϵ^(-1/n)
    symcurve = [symbolfunc(efac*cis(θ)) for θ in θv]
    plot(real(symcurve),imag(symcurve),"g")
    =#

    p
end

end
