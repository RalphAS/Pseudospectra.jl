# Caveat
This part of the documentation is a sketch. Please refer to examples
and test code.

# Usage

Typical use of the REPL interface is as follows:

```julia
using Pseudospectra
using Plots
A = your_matrix_generating_function()
setpsplotter()
gs = setgs()
ps_data = new_matrix(A)
options = Dict{Symbol,Any}()
driver!(ps_data,options,gs)
# modify `options` to concentrate on a region of interest
driver!(ps_data,options,gs)
```

This should show contour plots of ``\log_{10}(\epsilon)`` in the
vicinity of the spectrum - the standard display of a spectral
portrait. Eigenvalues of `A` are displayed as black points, if they
are available.

The `new_matrix` function constructs a stateful object with information
about `A` conducive to pseudospectral analysis; `driver!` then manages
the appropriate computations and plotting.


# Requirements

The integrated plotting capabilities require that the `Plots` and/or
`PyPlot` packages be installed. These are not formal package
requirements because much of the `Pseudospectra` package is useful
without plotting.

# Use without the integrated plotters

One can use the [`psa_compute`](@ref) function to "simply" evaluate
resolvent norms on a grid, as demonstrated in the `plots_minimal.jl`
script in the examples folder.
