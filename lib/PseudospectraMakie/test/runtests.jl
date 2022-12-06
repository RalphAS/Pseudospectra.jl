using Pseudospectra, PseudospectraMakie, Test, Pkg

parentdir = Pkg.pkgdir(Pseudospectra)

if length(ARGS) > 0
    tests = ARGS
else
    include(joinpath(parentdir, "test", "plot_tests.jl"))
end

psplotter = :Makie

using CairoMakie
@info("Plotting package is $psplotter")
setpsplotter(psplotter)
headless = isdefined(Main,:displayplots) ? !Main.displayplots : true
gs = setgs(headless=headless, savefigs=headless, fig_id=1)

for f in tests
    println("Testing $f functionality")
    fnam = f*".jl"
    include(joinpath(parentdir, "test", fnam))
end
