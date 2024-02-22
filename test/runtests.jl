using Pseudospectra, Test, LinearAlgebra
using Aqua

Aqua.test_all(Pseudospectra)

@testset "psa_compute" begin
    n = 32
    A = triu(rand(n, n))
    ax = [-5,5,-5,5]
    res = Pseudospectra.psa_compute(A, 40, ax, [], Dict{Symbol,Any}())
end

# more purely computational tests will go here...

# it would be grand if this worked:

# using CairoMakie
# include("../ext/PseudospectraMakie/test/runtests.jl")
