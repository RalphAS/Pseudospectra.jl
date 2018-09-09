#=
plot package initialization for Pseudospectra.jl

This file is part of Pseudospectra.jl.

SPDX-License-Identifier: BSD-3-Clause
License-Filename: LICENSES/BSD-3-Clause_Eigtool
=#
export setpsplotter, setgs, isheadless

const myplotter = Ref{Symbol}(:undef)

"""
    setpsplotter(plotter::Symbol=:default)

Select a plotting package for use with Pseudospectra.

Currently `:Plots` and `:PyPlot` are implemented.
Defaults to `:Plots` unless PyPlot is already imported without Plots.
"""
function setpsplotter(plotter::Symbol=:default)
    if plotter == :default
        if isdefined(Main, :PyPlot) && !isdefined(Main, :Plots)
            plotter = :PyPlot
        else
            plotter = :Plots
        end
    end
    if plotter == :Plots
        need2load = !isdefined(Main, :PseudospectraPlots)
        fnam = joinpath(@__DIR__,"PseudospectraPlots.jl")
    elseif plotter == :PyPlot
        need2load = !isdefined(Main, :PseudospectraMPL)
        fnam = joinpath(@__DIR__,"PseudospectraMPL.jl")
    else
        throw(ArgumentError("plotter argument must be :Plots or :PyPlot"))
    end
    if need2load
        Core.eval(Main, :(include($fnam)))
    end
    myplotter[]=plotter
    nothing
end

"""
    setgs(; headless=false, savefigs=true) => gs

Construct a `GUIState` for subsequent use by Pseudospectra functions.

Assumes plotting package has been chosen via `setpsplotter()`.
"""
function setgs(; headless=false, savefigs=true, fig_id=0)
    if myplotter[] ∈ [:default, :Plots]
        gs = Main.PseudospectraPlots.PlotsGUIState(nothing,fig_id,nothing,
                                                   headless=headless,
                                                   savefigs=savefigs)
    elseif myplotter[] == :PyPlot
        gs = Main.PseudospectraMPL.MPLGUIState(nothing,fig_id,nothing,
                                               headless=headless,
                                               savefigs=savefigs)
    else
        throw(ErrorException("use `setpsplotter` to establish plotting "
                             * "package for Pseudospectra first"))
    end
    return gs
end
