#=
Toeplitz matrix examples from Trefethen & Embree, Spectra & Pseudospectra
mainly ch. 7

include() this file, then invoke M.runme() with args of interest
=#
module M
using Pseudospectra, ToeplitzMatrices, Polynomials
using PyPlot


function runme(casename="limaçon",n=64,niter=100)

    # Puzzle: why does default work from Main/REPL, but not in a module?
    setpsplotter(:PyPlot)
    # setpsplotter(:Plots)
    # if using Plots, you should probably(?) set a backend first

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
        scale!(pert,ϵ/pn)
        A1 = A + pert
        ews = eigvals(A1)
        plot(real(ews),imag(ews),"r.")
#        scatter!(real(ews),imag(ews),markersize=2,c=:red,leg=false,
#        markerstrokealpha=0.0)
    end
    θv = linspace(0,2π,201)

    symcurve = [symbolfunc(cis(θ)) for θ in θv]
    plot(real(symcurve),imag(symcurve),"k")
    #    plot!(real(symcurve),imag(symcurve))

    #=
    # Sometimes perturbed ews are inside this curve.
    # Theory seems confused.

    efac = ϵ^(-1/n)
    symcurve = [symbolfunc(efac*cis(θ)) for θ in θv]
    plot(real(symcurve),imag(symcurve),"g")
    =#

    nothing
end

end
