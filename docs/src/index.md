```@meta
CurrentModule = NLS_Solver
```

```@setup session
using NLS_Solver
```

# Introduction

**NLS_Solver.jl** is a pure Julia package to solve unconstrained and
bound constrained nonlinear least squares problems.

- classical nonlinear least squares:
```math
\min\limits_{\theta} \frac{1}{2}\|r(\theta)\|^2
```
- bound constrained nonlinear least squares
```math
\min\limits_{\theta_l\le\theta\le\theta_u} \frac{1}{2}\|r(\theta)\|^2
```

## Usage

This package is easy to use. Problems are solved by calling a generic
`solve()` function of the form:

```julia
solve(problem, algorithm_conf)::algorithm_result
```

- **problem** is problem dependant data
- **algorithm_conf** defines the chosen algorithm
- **algorithm_result** is a structure that contains the
  result.

Detailed examples are provided belows:

```@contents
Pages = [
    "unconstrained_nls.md",
    "bound_constrained_nls.md",
    "nonlinear_regressions.md",
]
Depth = 3
```

## Detailed documentation

```@contents
Pages = [
    "api.md",
]
Depth = 3
```

## Algorithms and references 

For the moment two methods are implemented. The classical
Levenberg-Marquardt method for unconstrained problems and a modified
Levenberg-Marquardt method for bound constrained problems.

These implementations are mainly based on these references:

- Madsen, N. (). Methods For Non-Linear Least Squares Problems.
  [imm3215.pdf](http://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf) 
  
  This reference is a must read about nonlinear least squares
  algorithms.

- Nielsen, H. B., & others, (1999). Damping parameter in marquardt's
  method. [tr05_99.pdf](http://www2.imm.dtu.dk/documents/ftp/tr99/tr05_99.pdf)

  This second reference provides some details about damping
  parameter. It also gives an useful list of test functions
  (unconstrained case).

- Kunisch, K., & Rendl, F. (2003). An infeasible active set method for
  quadratic problems with simple bounds. SIAM Journal on Optimization,
  14(1), 35–52. [epubs.siam.org](http://dx.doi.org/10.1137/s1052623400376135)

  This reference presents an efficient method to solve bound
  constrained quadratic problems.

