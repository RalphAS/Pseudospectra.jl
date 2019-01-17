#=
plot package initialization for Pseudospectra.jl

This file is part of Pseudospectra.jl.

SPDX-License-Identifier: BSD-3-Clause
License-Filename: LICENSES/BSD-3-Clause_Eigtool
=#
export setpsplotter, getpsplotter, setgs, isheadless, defaultgs

const _currentplotter = Ref{Symbol}(:undef)

const _available_plotters = Symbol[:PyPlot, :Plots]
const _enabled_plotters = Dict{Symbol,Bool}(:PyPlot => false, :Plots => false)

const _currentgs = Ref{Union{Nothing,GUIState}}(nothing)

function link_pyplot()
    include("PseudospectraMPL.jl")
    _enabled_plotters[:PyPlot] = true
end

function link_plots()
    include("PseudospectraPlots.jl")
    _enabled_plotters[:Plots] = true
end

getpsplotter() = _currentplotter[]

"""
    setpsplotter(plotter::Symbol=:default)

Select a plotting package for use with Pseudospectra.

Currently `:Plots` and `:PyPlot` are implemented.
Defaults to `:Plots` unless PyPlot is already imported without Plots.
"""
function setpsplotter(plotter::Symbol=:default)
    if plotter == :default
        if _enabled_plotters[:PyPlot] && !_enabled_plotters[:Plots]
            plotter = :PyPlot
        else
            plotter = :Plots
        end
    end
    if plotter ∉ _available_plotters
        throw(ArgumentError("plotter argument must be :Plots or :PyPlot"))
    end
    if !_enabled_plotters[plotter]
        error("selected or default plotter '$plotter' is not enabled")
    end
    _currentplotter[]=plotter
    nothing
end

"""
    setgs(; headless=false, savefigs=true) => gs

Construct a `GUIState` for subsequent use by Pseudospectra functions.

Assumes plotting package has been chosen via `setpsplotter()`.
"""
function setgs(; headless=false, savefigs=true, fig_id=0)
    if _currentplotter[] ∈ [:default, :Plots]
        gs = PseudospectraPlots.PlotsGUIState(nothing,fig_id,nothing;
                                                   headless=headless,
                                                   savefigs=savefigs)
    elseif _currentplotter[] == :PyPlot
        gs = PseudospectraMPL.MPLGUIState(nothing,fig_id,nothing;
                                               headless=headless,
                                               savefigs=savefigs)
    else
        throw(ErrorException("use `setpsplotter` to establish plotting "
                             * "package for Pseudospectra first"))
    end
    _currentgs[] = gs
    return gs
end

"""
provides a default value for the `GUIState` object used by the Pseudospectra
package.  Attempts to initialize such an object if necessary.
"""
function defaultgs()
    if _currentgs[] === nothing
        if _currentplotter[] == :undef
            setpsplotter()
        end
        _currentgs[] = setgs()
    end
    _currentgs[]
end
