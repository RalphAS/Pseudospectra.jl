# Pseudospectra

[![Build Status](https://travis-ci.org/RalphAS/Pseudospectra.jl.svg?branch=master)](https://travis-ci.org/RalphAS/Pseudospectra.jl)
[![Coverage Status](http://codecov.io/github/RalphAS/Pseudospectra.jl/coverage.svg?branch=master)](http://codecov.io/github/RalphAS/Pseudospectra.jl?branch=master)

# Introduction
Pseudospectra is a Julia package for computing pseudospectra of
non-symmetric matrices, and plotting them along with eigenvalues
("spectral portraits"). Some related computations and plots are
also provided.

## Mathematical background
The ϵ-pseudospectrum of a matrix `A` is the locus (in the complex plane) of
eigenvalues of perturbations `A+E` where `norm(E) < ϵ`. Specifically,
this package works with the unweighted 2-norm.

Among other things, pseudospectra:
* elucidate transient behavior hidden to eigen-analysis, and
* indicate the utility of eigenvalues extracted via iterative methods like `eigs`.

See [the Pseudospectra gateway](http://www.cs.ox.ac.uk/pseudospectra/intro.html)
for details, references, and more.

## Package context
Pseudospectra (along with the QML-based GUI, in the forthcoming PseudospectraView
package) is essentially a translation of the acclaimed MATLAB-based EigTool
([homepage here](http://www.comlab.ox.ac.uk/pseudospectra/eigtool)),
code now hosted [on GitHub](https://github.com/eigtool/eigtool).

No endorsement or promotion of Pseudospectra.jl by the authors of EigTool
is implied.

Specific documentation for Pseudospectra is a work in progress. See the
examples and tests.

## Note on packaging/requirements
The plotting interface is somewhat schizophrenic. Drivers are included
for Plots.jl and/or PyPlot.jl (i.e., PyPlot is a useful back end for
Plots as used here; other Plots backends have been partially tested).

Although this package is designed with an eye to plotting results,
the computational routines are usable without a plotting package,
To avoid forcing potential users to install a particular one, none are
specified in the formal package requirements.  The `setgs` function
can import one conditionally.

Some functions used for examples require other packages. They should
give a useful complaint if invoked without that support.

# Installation
Until this project is mature enough to be registered, install by

```julia
Pkg.clone("https://github.com/RalphAS/Pseudospectra.jl")
```

# Basics
Minimal use of the REPL interface is as follows:

```julia
using Pseudospectra
A = your_matrix_generating_function()
setpsplotter()
gs = setgs()
ps_data = new_matrix(A)
options = Dict{Symbol,Any}()
driver!(ps_data,options,gs)
```

More realistic uses are (as always) in the examples and test folders.

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
The source files are under the BSD license, in accordance with derivation
from EigTool, except for `xeigs.jl` which derives from Julialang code and
is under an MIT license.
