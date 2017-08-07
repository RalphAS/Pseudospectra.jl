module M
# numerical range
# This test is based on T&E chap. 43
using Pseudospectra
using Base.Test

include("testplotter.jl")
η = 0.015
A = Pseudospectra.advdiff(100,η)
opts = Dict{Symbol,Any}(:showfov => true, :fov_npts => 200,
                        :ax => [-60,20,-50,50])

ps_data = new_matrix(A,opts)
driver!(ps_data,opts,gs)

@test isapprox(maximum(real(ps_data.ps_dict[:fov])),-η*π^2,atol=1e-6)
end
