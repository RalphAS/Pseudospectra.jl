#=
plot package initialization for Pseudospectra.jl

This file is part of Pseudospectra.jl.

SPDX-License-Identifier: BSD-3-Clause
License-Filename: LICENSES/BSD-3-Clause_Eigtool
=#
export setpsplotter, getpsplotter, setgs, isheadless, defaultgs

const _currentplotter = Ref{Symbol}(:undef)

const _available_plotters = Symbol[:PyPlot, :Plots,
                                   :Makie,
                                   ]
const _enabled_plotters = Dict{Symbol,Bool}()
const _guistates = Dict{Symbol,Any}()

const _currentgs = Ref{Union{Nothing,GUIState}}(nothing)

# We could restrict args, but footguns are the Julian way
function _register_plotter(s::Symbol, guistatefn)
    _enabled_plotters[s] = true
    _guistates[s] = guistatefn
end

getpsplotter() = _currentplotter[]

"""
    setpsplotter(plotter::Symbol)

Select a plotting package for use with Pseudospectra.

This is typically invoked automatically upon loading one of the associated packages.
"""
function setpsplotter(plotter::Symbol)
    if ! get(_enabled_plotters,plotter,false)
        throw(ArgumentError("Selected plotter '$plotter' is not enabled;\nan appropriate package (e.g. Pseudospectra$(plotter)) must be loaded first."))
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
    fn = get(_guistates, _currentplotter[], nothing)
    if fn === nothing
        throw(ErrorException("use `setpsplotter` to establish plotting "
                             * "package for Pseudospectra first"))
    end
    gs = fn(nothing,fig_id,nothing; headless=headless, savefigs=savefigs)
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
