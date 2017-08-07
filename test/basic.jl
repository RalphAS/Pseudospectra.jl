module M
using Pseudospectra
# using Base.Test

include("testplotter.jl")

A = Pseudospectra.grcar(40)

opts = Dict{Symbol,Any}()
ps_data = new_matrix(A,opts)
driver!(ps_data,opts,gs)
modeplot(ps_data,gs,0,0.5+2.0im)
modeplot(ps_data,gs,1,0.5+2.0im)

# @test res[1] == exactval1
# @test res[2] ≈ approxval2

end
