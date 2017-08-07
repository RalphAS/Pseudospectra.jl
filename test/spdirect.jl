module M
# 400x400 real sparse direct
using Pseudospectra
# using Base.Test
include("testplotter.jl")

# for sparse direct, must provide axes
opts = Dict{Symbol,Any}(:ax => [-0.5,1.0,1.5,3.0],:npts=>60,
                        :direct => true, :sparse_direct => true)
A = sparse(Pseudospectra.grcar(400))

ps_data = new_matrix(A,opts)
driver!(ps_data,opts,gs)

end
