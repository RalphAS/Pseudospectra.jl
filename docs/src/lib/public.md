# Public interface

## Importing a matrix or operator

```@docs
new_matrix
```

## Setting up graphics subsystem

```@docs
setpsplotter

setgs
```

## Driver function
After a matrix has been ingested into a `PSAStruct` and the graphics
subsystem has been established, the following function will compute
pseudospectra and plot a spectral portrait:

```@docs
driver!
```

## Pseudospectra computation

```@docs
psa_compute
```

## Eigen/Pseudo-mode computation and plotting

```@docs
modeplot
```

## Other computations

```@docs
psa_radius

psa_abscissa

numerical_range

numerical_abscissa
```

## Other plots

```@docs
Pseudospectra.surfplot

mtxexpsplot

mtxpowersplot
```
