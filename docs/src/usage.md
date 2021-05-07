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
resolvent norms on a grid, as follows:
```julia
using LinearAlgebra, Pseudospectra
# matrix must be complex to get a true Schur decomposition
A = your_matrix_generator() .+ 0im
Tschur, U, eigA = schur(A)
ax = [-3,3,-3,3] # domain of results
npts = 100  # mesh size
opts = Dict{Symbox,Any}()
Z, x, y, levels = psa_compute(Tschur, npts, ax, eigA, opts)
```

## Use with extended precision

The dense-matrix methods should work with matrix element types of extended precision
such as `Float128` and `BigFloat`, but they require linear algebra methods which
are implemented in the `GenericLinearAlgebra` and `GenericSchur` packages.

```julia
using Pseudospectra, Plots, GenericLinearAlgebra, GenericSchur, Quadmath
Aq = Float128.(A_matrix)
spectralportrait(A)
```

Note that linear algebra with non-BLAS types is expensive, so computation takes many
minutes even for moderate-sized matrices.

## Parallelism

The computation of spectral portraits lends itself to parallel processing.
Multiprocessing is only implemented for some algorithm variants.

The core algorithm can be set to distribute the `Z`-grid over Julia
threads with the `threaded` entry in `options`.  This is especially useful for
extended precision, and for working on fine grids with BLAS
types. (One may need to adjust thread usage in the BLAS to balance if
the library is aggressive).

Unfortunately `BigFloat` operations currently make heavy use of the heap, invoking
the thread-unfriendly garbage collector, so for extended precision `Float128` may
be more practical.

For larger matrices, parallelism may best be left to the underlying linear algebra libraries.
