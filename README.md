# Pseudospectra

[![GitHub CI Build Status](https://github.com/RalphAS/Pseudospectra.jl/workflows/CI/badge.svg)](https://github.com/RalphAS/Pseudospectra.jl/actions)
[![Coverage Status](http://codecov.io/github/RalphAS/Pseudospectra.jl/coverage.svg?branch=master)](http://codecov.io/github/RalphAS/Pseudospectra.jl?branch=master)
[![Documentation/dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://RalphAS.github.io/Pseudospectra.jl/dev)
[![Documentation/dev](https://img.shields.io/badge/docs-stable-blue.svg)](https://RalphAS.github.io/Pseudospectra.jl/stable)

# Introduction
Pseudospectra is a Julia package for computing pseudospectra of
non-symmetric matrices, and plotting them along with eigenvalues
("spectral portraits"). Some related computations and plots are
also provided.

## Mathematical background
Whereas the spectrum of a matrix is the set of its eigenvalues,
a pseudospectrum is the set of complex numbers "close" to the spectrum
in some practical sense.

More precisely, the ϵ-pseudospectrum of a matrix `A`, `σ_ϵ(A)`, is the set of
complex numbers `λ` such that
* `λ` is an eigenvalue of some matrix `A+E`, where the norm of the perturbation `‖E‖ < ϵ`, or
* the resolvent norm `‖(A-λI)^(-1)‖ > 1/ϵ`,

(the definitions are equivalent). This sense of "closeness" is trivial for
Hermitian matrices, but interesting for others.
Specifically, this package is currently limited to the unweighted 2-norm.

Among other things, pseudospectra:
* elucidate transient behavior hidden to eigen-analysis, and
* indicate the utility of eigenvalues extracted via iterative methods like `eigs` (from the Arpack package).

See [the Pseudospectra gateway](http://www.cs.ox.ac.uk/pseudospectra/intro.html)
for details, references, and more.

## Aside: the simple interface
To study a moderate-sized matrix with minimal user effort,
follow this example:
```julia
using Plots, Pseudospectra, LinearAlgebra
n=150
B=diagm(1 => fill(2im,n-1), 2 => fill(-1,n-2), 3 => fill(2,n-3), -2 => fill(-4,n-2), -3 => fill(-2im, n-3))
spectralportrait(B)
```

![example figure](https://user-images.githubusercontent.com/18298838/55284298-c4213100-5341-11e9-8718-514acdf3ab9e.png)

The figure shows a section of the complex plane with eigenvalues and contours
of `log10(ϵ)`.

## Package context
Pseudospectra.jl (along with associated graphical packages) is largely a translation of
the acclaimed MATLAB-based EigTool
([homepage here](http://www.comlab.ox.ac.uk/pseudospectra/eigtool)),
code now hosted [on GitHub](https://github.com/eigtool/eigtool).

No endorsement or promotion of Pseudospectra.jl by the authors of EigTool
is implied.

Specific documentation for Pseudospectra is a work in progress; a draft is
available [here](https://RalphAS.github.io/Pseudospectra.jl/stable). See the
examples and tests for more.


# Installation
This package is included in the General registry,
so the normal `Pkg` commands to add `Pseudospectra` suffice.

## Note on packaging/requirements
Although this package is designed with an eye to plotting results,
the computational routines are usable without a plotting package,

Drivers are included for Plots.jl, Makie.jl, and PyPlot.jl. The drivers will be loaded
automatically if the corresponding plotting package(s) is/are loaded (via Requires or as
package extensions).

Some functions used for examples require other packages. They should
give a useful complaint if invoked without that support.

# Basic usage
Minimal use of the REPL interface is as follows:

```julia
using Plots
using Pseudospectra
A = your_matrix_generating_function()
ps_data = new_matrix(A)
driver!(ps_data)
# modify, e.g., for higher resolution
options = Dict{Symbol,Any}(:npts => 100)
driver!(ps_data,options)
```

This should show a contour plot of `log10(ϵ)` in the vicinity of the spectrum,
which is the standard display of a spectral portrait.
More elaborate capabilities are exhibited (as always) in the examples and
test folders.

# Disclaimer
This software is provided by the copyright holders and contributors "as is" and
any express or implied warranties, including, but not limited to, the implied
warranties of merchantability and fitness for a particular purpose are
disclaimed. In no event shall the copyright holder or contributors be liable for
any direct, indirect, incidental, special, exemplary, or consequential damages
(including, but not limited to, procurement of substitute goods or services;
loss of use, data, or profits; or business interruption) however caused and
on any theory of liability, whether in contract, strict liability, or tort
(including negligence or otherwise) arising in any way out of the use of this
software, even if advised of the possibility of such damage.

# Note on licensing
Most of the package is under a BSD license, in accordance with derivation
from EigTool. See individual source files for exceptions.
