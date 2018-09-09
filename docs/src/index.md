# Pseudospectra.jl Documentation

Pseudospectra is a Julia package for computing pseudospectra of
non-symmetric matrices, and plotting them along with eigenvalues
("spectral portraits"). Some related computations and plots are
also provided.

## Introduction
Whereas the spectrum of a matrix is the set of its eigenvalues,
a pseudospectrum is the set of complex numbers "close" to the spectrum
in some practical sense.

More precisely, the ϵ-pseudospectrum of a matrix `A`, ``\sigma_{\epsilon}(A)``,
is the set of complex numbers ``\lambda`` such that
* ``\lambda`` is an eigenvalue of some matrix ``A+E``, where the perturbation ``E`` is small: ``\|E\| < \epsilon``
* the resolvent at ``\lambda`` has a large norm: ``\|(A-λI)^{-1}\| > 1/\epsilon``,

(the definitions are equivalent).
Specifically, this package is currently limited to the unweighted 2-norm.

Among other things, pseudospectra:
* elucidate transient behavior hidden to eigen-analysis, and
* indicate the utility of eigenvalues extracted via iterative methods like `eigs`.

This package facilitates computation, display, and investigation of
the pseudospectra of matrices and some other representations of linear
operators.

## Spectral portraits
It is customary to display pseudospectra as contour plots of the logarithm
of the inverse of the resolvent norm ``\epsilon = 1/\|(A-zI)^{-1}\|`` for ``z``
in a subset of the complex plane. Thus ``\sigma_{\epsilon}(A)`` is the union
of the interiors of such contours. Such plots, sometimes called
*spectral portraits*, are the most prominent product of this package.

## Credit
Pseudospectra.jl is essentially a translation of the acclaimed MATLAB-based EigTool
([homepage here](http://www.comlab.ox.ac.uk/pseudospectra/eigtool))

## References
* [The Pseudospectra gateway](http://www.cs.ox.ac.uk/pseudospectra/intro.html).
* L.N. Trefethen and M.Embree, *Spectra and Pseudospectra; The Behavior of Nonnormal Matrices and Operators*, Princeton 2005,
