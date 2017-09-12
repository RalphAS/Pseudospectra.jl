# Caveat
This part of the documentation is a placeholder. Please refer to examples
and test code.

# Usage

Typical use of the REPL interface is as follows:

```julia
using Pseudospectra
A = your_matrix_generating_function()
setpsplotter()
gs = setgs()
ps_data = new_matrix(A)
options = Dict{Symbol,Any}()
driver!(ps_data,options,gs)
```

This should show a contour plot of ``\log_{10}(\epsilon)`` in the vicinity of the spectrum,
which is the standard display of a spectral portrait. Eigenvalues of `A` are
displayed as black points.

# Requirements
The plotting capabilities require that the `Plots` and/or `PyPlot` packages
be installed. These are not formal package requirements because much of
the `Pseudospectra` package is useful without plotting.
