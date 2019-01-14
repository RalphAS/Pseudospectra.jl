# minimal use of the package, user manages all the plotting
# PyPlot version

module M
using Pseudospectra
using PyPlot

A = randn(100,100)
# we want a true Schur decomposition, so force complexity
Tschur,U,eigA  = schur(A .+ 0im)
ax = [-12,12,-12,12]
npts = 100
opts = Dict{Symbol,Any}(:real_matrix => true)
Z,x,y,levels,err,Tproj,eigAproj,algo = psa_compute(Tschur,npts,
                                                   ax,eigA,opts)
plot(real(eigA),imag(eigA),"k.")
if isempty(levels)
    contour(x,y,log10.(Z),extend="both",linewidths=2)
else
    contour(x,y,log10.(Z),levels=levels,extend="both",linewidths=2)
end
colorbar()
end
