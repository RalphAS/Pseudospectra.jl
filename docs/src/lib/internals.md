# Internals

## Data structures

```@docs
Portrait

PSAStruct
```

## Arnoldi iteration

```@docs
Pseudospectra.xeigs

ArpackOptions
```

## Algorithm selection
IRAM means implicitly restarted Arnoldi method, using the ARPACK implementation.
"Large" and "small" are determined by constants which might be made
configurable someday.

### Dense square matrix
* If large and not `opts[:direct]`, project via IRAM and go to rectangular
  branch.
* If small, use SVD.
* Otherwise, use inverse Lanczos.
  * If a Schur decomposition is available, use that to simplify the solvers.

### Sparse square matrix
* If small, and not `opts[:keep_sparse]`, convert to dense then use above logic.
* If `opts[:direct]`, use inverse Lanczos. This involves a sparse factorization
  for every grid point, so is slow.
* Otherwise project via IRAM and go to rectangular branch.

### Dense rectangular matrix
* If not Hessenberg, factor first.
  * If very tall, factor by QR, otherwise QZ (generalized Schur).
* Use QR for resolvents, then inverse Lanczos.
  * Use a simplified QR for large Hessenberg.
