using Pseudospectra, Test, LinearAlgebra
using Aqua

include("testplotter.jl")

Aqua.test_all(Pseudospectra)

@testset "psa_compute" begin
    n = 32
    A = triu(rand(n, n))
    ax = [-5,5,-5,5]
    res = Pseudospectra.psa_compute(A, 40, ax, [], Dict{Symbol,Any}())
end

# more purely computational tests will go here...

if psplotter == :None
    exit()
end

# define the "tests" array:
include("plot_tests.jl")

for f in tests
    println("Testing $f functionality")
    fnam = f*".jl"
    include(fnam)
end
