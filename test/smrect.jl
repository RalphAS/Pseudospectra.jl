module M
# 40x39 small rectangular
using Pseudospectra
# using Base.Test
include("testplotter.jl")

G = Pseudospectra.grcar(40)
A = G[:,1:end-1]
# for a rectangular mtx, must provide axes
opts = Dict{Symbol,Any}(:ax => [-1.0,3.0,-3.0,3.0],:npts=>60)

ps_data = new_matrix(A,opts)
driver!(ps_data,opts,gs)

end
