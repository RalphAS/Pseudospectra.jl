# Pseudospectra.jl Documentation

Pseudospectra is a Julia package for computing pseudospectra of
non-symmetric matrices, and plotting them along with eigenvalues
("spectral portraits"). Some related computations and plots are
also provided.

## Introduction
The ϵ-pseudospectrum of a matrix `A`, ``\sigma_{\epsilon}(A)``, is the set of complex
numbers ``\lambda`` such that
* there is a perturbation ``E`` where ``\|E\| < \epsilon`` such that ``\lambda`` is an eigenvalue of the matrix ``A+E``,
* the resolvent norm ``\|(A-λI)^{-1}\| > 1/\epsilon``,

(the definitions are equivalent).
Specifically, this package is currently limited to the unweighted 2-norm.

Among other things, pseudospectra:
* elucidate transient behavior hidden to eigen-analysis, and
* indicate the utility of eigenvalues extracted via iterative methods like `eigs`.
