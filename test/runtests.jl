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

# define the "tests" array:
include("plot_tests.jl")

using CairoMakie
psplotter = :Makie
@info("Plotting package is $psplotter")
setpsplotter(psplotter)
headless = isdefined(Main,:displayplots) ? !Main.displayplots : true
gs = setgs(headless=headless, savefigs=headless, fig_id=1)

for f in tests
    println("Testing $f functionality")
    fnam = f*".jl"
    include(fnam)
end
