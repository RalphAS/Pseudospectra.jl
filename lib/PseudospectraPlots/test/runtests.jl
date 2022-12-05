using Pseudospectra, PseudospectraPlots, Test, Pkg

tests = ["basic","threads","arpack1","medcplx","medrect","smrect","spdirect",
          "projection_and_transient","radius_abscissa","negproj","eigsnoconv",
         "numrange", "linmap", "big", "power_transient", "swmethod", "zoom",
          ]

if length(ARGS) > 0
    tests = ARGS
end

psplotter = :Plots

using Plots
@info("Plotting package is $psplotter")
setpsplotter(psplotter)
headless = isdefined(Main,:displayplots) ? !Main.displayplots : true
gs = setgs(headless=headless, savefigs=headless, fig_id=1)

for f in tests
    println("Testing $f functionality")
    fnam = f*".jl"
    include(joinpath(Pkg.pkgdir(Pseudospectra),"test",fnam))
end
