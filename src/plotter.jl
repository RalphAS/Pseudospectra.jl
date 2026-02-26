#=
plot package initialization for Pseudospectra.jl

This file is part of Pseudospectra.jl.

SPDX-License-Identifier: BSD-3-Clause
License-Filename: LICENSES/BSD-3-Clause_Eigtool
=#
export defaultgs, getpsplotter, isheadless, setgs, setpsplotter

const _currentplotter = Ref{Symbol}(:undef)

const _available_plotters = Symbol[:PyPlot, :PythonPlot, :Plots,
                                   :Makie,
                                   ]
const _enabled_plotters = Dict{Symbol, PSAPlotter}()
const _guistates = Dict{Symbol,Any}()

const _currentgs = Ref{Union{Nothing,GUIState}}(nothing)

# We could restrict args, but footguns are the Julian way
function _register_plotter(s::Symbol, guistatefn, pp)
    _enabled_plotters[s] = pp
    _guistates[s] = guistatefn
end

function getpsplotter()
    if _currentplotter[] == :undef
        ks = keys(_enabled_plotters)
        nk = length(ks)
        if nk == 0
            throw(ErrorException("No plotting interface is enabled.\n"
                                 * "An appropriate package (i.e. Plots, PythonPlot,"
                                 * "or Makie) must be loaded."))
        else
            # we can't index Dict keys, but we can iterate
            for k in ks
                setpsplotter(k)
                if nk > 1
                    @info "Using $k as plotting package for Pseudospectra"
                end
                break
            end
        end
    end
    return _enabled_plotters[_currentplotter[]]
end

"""
    setpsplotter(plotter::Symbol)

Select a plotting package for use with Pseudospectra.

This is typically invoked automatically upon loading one of the associated packages,
and may also be used to select one if several are loaded.
"""
function setpsplotter(plotter::Symbol)
    if get(_enabled_plotters, plotter, nothing) === nothing
        throw(ArgumentError("Selected plotter '$plotter' is not enabled;\n"
                            * "an appropriate package (i.e. Plots, PythonPlot, or Makie)"
                            * "must be loaded first."))
    end
    _currentplotter[]=plotter
    nothing
end

"""
    setgs(; headless=false, savefigs=true) => gs

Construct a `GUIState` for subsequent use by Pseudospectra functions.

Assumes plotting package has been chosen via `setpsplotter()`.
"""
function setgs(; drawcmd=nothing, headless=false, savefigs=true, fig_id=0, kwargs...)
    fn = get(_guistates, _currentplotter[], nothing)
    if fn === nothing
        throw(ErrorException("use `setpsplotter` to establish plotting "
                             * "package for Pseudospectra first"))
    end
    gs = fn(nothing, fig_id, drawcmd; headless=headless, savefigs=savefigs, kwargs...)
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
            getpsplotter()
        end
        _currentgs[] = setgs()
    end
    _currentgs[]
end
