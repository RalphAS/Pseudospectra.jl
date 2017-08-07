module M
# Test for graceful return if xeigs does not converge
# This uses the tridiagonal Toeplitz matrix from T&E chap. 3,
# which is nasty enough to fail Arnoldi.

# CHECKME:
# The start vector for xeigs is pseudorandom.
# It is conceivable that some realizations might cause this test to fail.

using Pseudospectra
using Base.Test

include("testplotter.jl")

opts = Dict{Symbol,Any}(:npts=>100,:levels=>-12:1:-1,:ax=>[-1.5,1.5,-1,1],
                        :direct=>false)

N = 128
A = diagm(ones(N-1),1) + (1/4)* diagm(ones(N-1),-1)

ps_data = new_matrix(A,opts)
driver!(ps_data,opts,gs)

@test !ps_data.zoom_list[1].computed

end
