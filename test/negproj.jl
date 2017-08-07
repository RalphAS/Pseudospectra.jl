module M
# complex dense, projected away from origin
using Pseudospectra
# using Base.Test
include("testplotter.jl")

# for a rectangular mtx, must provide axes
opts = Dict{Symbol,Any}(:ax => [-1.2,1.2,-1.1,1.1], :npts => 60,
                        :levels => -5:0.5:-1, :proj_lev => -5.0)
A = Pseudospectra.landau_fox_li(225,16)

ps_data = new_matrix(A,opts)
driver!(ps_data,opts,gs)

end
