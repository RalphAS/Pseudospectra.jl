# Usage

To use the plotting capabilities, one should first explicitly import
either `Plots` or `PyPlot`. (An interface to `Makie` is in progress.)

A simple interface for modestly-sized matrices is available through the
[`spectralportrait`](@ref) function.

Typical "advanced" use of the package is as follows:

```julia
using Pseudospectra
using Plots
A = your_matrix_generating_function()
ps_data = new_matrix(A)
driver!(ps_data)
options = Dict{Symbol,Any}()
# modify `options` to concentrate on a region of interest, increase resolution, etc.
driver!(ps_data,options)
```

This should show contour plots of ``\log_{10}(\epsilon)`` in the
vicinity of the spectrum - the standard display of a spectral
portrait. Eigenvalues of `A` are displayed as black points, if they
are available.

The `new_matrix` function constructs a stateful object with information
about `A` conducive to pseudospectral analysis; `driver!` then manages
the appropriate computations and plotting.


## Requirements

The integrated plotting capabilities require that the `Plots` and/or
`PyPlot` packages be installed. These are not formal package
requirements because much of the `Pseudospectra` package is useful
without plotting.

## Annoyances

(These are problems with other packages that may arise here.  If you
find problems with Pseudospectra.jl itself, please file an issue.)

If you are using IJulia (viz. Jupyter) and Plots, you may need to
issue `inline()` to get the plots to appear.

Some of the supplementary plotting functions use LaTeX strings which
are displayed incorrectly by certain plotting backends.

## Use without the integrated plotters

One can use the [`psa_compute`](@ref) function to "simply" evaluate
resolvent norms on a grid, as demonstrated in the `plots_minimal.jl`
script in the examples folder.
