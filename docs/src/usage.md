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

The [`new_matrix`](@ref) function constructs a stateful object with information
about `A` conducive to pseudospectral analysis; [`driver!`](@ref) then manages
the appropriate computations and plotting.


## Requirements

The integrated plotting capabilities require that the `Plots` and/or
`PyPlot` packages be installed. These are not formal package
requirements because much of the `Pseudospectra` package is useful
without plotting.

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

The dense-matrix methods work with matrix element types of extended precision
such as `Float128` (from the `Quadmath` package) and `BigFloat`; they do require
linear algebra methods which
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
Multiprocessing is implemented here for some algorithm variants.

The core dense-matrix algorithms can be set to distribute the `Z`-grid over Julia
threads with the `threaded` entry in `options`.  This is especially useful for
extended precision, and for working on fine grids with BLAS
types. (One may need to adjust thread usage in the BLAS to balance if
the library is aggressive).

Unfortunately `BigFloat` operations currently make heavy use of the heap, invoking
the thread-unfriendly garbage collector, so for extended precision `Float128` may
be more practical.

## Large matrices

Currently, the methods for handling large matrices are those ported from Eigtool.
They can be invoked via options to the user-facing computational functions.

### Direct-sparse

This is a Lanczos scheme which does a sparse LU decomposition for each point in the `Z` grid.

### Regional subspace projection

This uses one Schur factorization of the matrix and is not always reliable. [More
documentation needed]

### Arnoldi projection

This uses a specialized wrapper of the ARPACK library, and may also be used for linear maps
not explicitly represented as matrices. [More documentation needed]
