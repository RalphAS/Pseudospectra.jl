# Prepare for plotting in test environment

# Default is to run offline tests, producing PNG files in temp. dir.,
# with default plotting package (Plots w/pyplot backend).

# To test other packages/backends, set them up first and define `psplotter`;
# e.g. do
# `using Plots; gr(); psplotter = :Plots`
# before running the test suite.

# To get figures on screen instead of PNG files, do
# `displayplots=true`
# before running the test suite.

if isdefined(Main,:psplotter)
    psplotter = Main.psplotter
    # warning: assume plotting package is properly initialized
else
    psplotter = Symbol(get(ENV,"PSPLOTTER","Plots"))
    if psplotter == :Plots
        using Plots
#        pyplot() # seems most reliable for png output
    elseif psplotter == :PyPlot
        using PyPlot
    elseif psplotter == :Makie
        using GLMakie
    else
        error("invalid ENV[\"PSPLOTTER\"]")
    end
end
@info("Plotting package is $psplotter")
setpsplotter(psplotter)
headless = isdefined(Main,:displayplots) ? !Main.displayplots : true
gs = setgs(headless=headless, savefigs=headless, fig_id=1)
