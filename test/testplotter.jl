# Prepare for plotting in test environment

# sets global variables: psplotter, gs

# This is a partial script rather than a function definition
# because mucking with the global state needs to be done at top level
# to avoid world age discrepancies.


"""
Set up the plotting backend for Pseudospectra in a testing run.

Default is to run offline tests, producing PNG files in temp. dir.,
with default plotting package (Plots w/default backend).


To test other packages/backends, either
a) set them up first and define `psplotter`;
  e.g. do
    `using Plots; gr(); psplotter = :Plots`
  before running the test suite,
or
b) use the envar PSPLOTTER and hope the environment cooperates.

To get figures on screen instead of PNG files, do
`displayplots=true` in the Main module
before running the test suite.
"""

const DEFAULT_PLOTTER = "Plots"

if isdefined(Main,:psplotter)
    psplotter = Main.psplotter
    # warning: assume plotting package is properly initialized
else
    let envstr = get(ENV, "PSPLOTTER", DEFAULT_PLOTTER)
        p = occursin("Makie", envstr) ? :Makie : Symbol(envstr)
        global psplotter = p
    end
    if psplotter == :Plots
        using Plots
        # pythonplot() # seems most reliable for png output but adds baggage
    elseif psplotter == :PythonPlot
        using PythonPlot # would be @eval'd if in a function
    elseif psplotter == :PyPlot
        using PyPlot
    elseif psplotter == :Makie
        using CairoMakie
    elseif psplotter != :None
        error("invalid ENV[\"PSPLOTTER\"];" *
              " expected Plots, PythonPlot, PyPlot, *Makie or None")
    end
end

if psplotter == :None
    @info("Skipping plotting functionality for this run.")
    gs = nothing
else
    @info("Plotting package is $psplotter.")
    setpsplotter(psplotter)
    headless = isdefined(Main,:displayplots) ? !Main.displayplots : true
    gs = setgs(headless=headless, savefigs=headless, fig_id=1)
end
