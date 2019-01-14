# minimal use of the package, user manages all the plotting

using Pseudospectra
using Plots

A = randn(100,100)
# we want a true Schur decomposition, so force complexity
Tschur,U,eigA  = schur(A .+ 0im)
ax = [-12,12,-12,12]
npts = 40
opts = Dict{Symbol,Any}(:real_matrix => true)
Z,x,y,levels,err,Tproj,eigAproj = psa_compute(Tschur,npts,
                                              ax,eigA,opts)
p = scatter(real(eigA),imag(eigA),color="black",label="")
if isempty(levels)
    contour!(p,x,y,log10.(Z),linewidth=2)
else
    contour!(p,x,y,log10.(Z),levels=levels,linewidth=2)
end
