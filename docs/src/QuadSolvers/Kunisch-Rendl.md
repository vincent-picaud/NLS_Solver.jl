# Kunisch-Rendl solver

## Introduction

This file implements the [Kunisch-Rendl method](https://epubs.siam.org/doi/10.1137/S1052623400376135). This methods is a primal-dual active set method for quadratic problems with bound constraints.

**Attention:** this method is not always convergent. But:
- in our context this is not a problem as we can enforce strong
  diagonal dominance of the matrix and the method is convergent.
- if required there is a modification of the algorithm,
  [Hungerlaender-Rendl
  paper](http://www.optimization-online.org/DB_HTML/2016/09/5644.html),
  that insures convergence.

## Documentation

```@autodocs
Modules = [NLS_Solver]
Pages = ["QuadSolvers/Kunisch-Rendl.jl"]
```
