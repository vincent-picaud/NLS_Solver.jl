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
    "nonlinear_regressions.md",
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

## More usage examples

This package is used by the
[NLS_Fit.jl](https://github.com/vincent-picaud/NLS_Fit.jl) package
which is dedicated to peak fitting in spectrometry.




To use this
method call [`LevenbergMarquardt_Conf`](@ref). The other implemented
method is a modification of the Levenberg-Marquardt where the inner
quadratic problem is solved by the Kunisch-Rendl method to handle
bound constraints. To use this method call [`LevenbergMarquardt_BC_Conf`](@ref). 
Please read this package
[README.org](https://github.com/vincent-picaud/NLS_Solver.jl) for
extra details.

The two following examples illustrate how to solve a classical
nonlinear least squares problem and a bound constrained one.

In both cases, the objective function is the Rosenbrock function:
```math
(1-θ_1)^2 + 100 (θ_2 - θ_1^2)^2
```
The classical problem is solved as follows:

```@example session
# define the objective methods
nls = create_NLS_problem_using_ForwardDiff(2 => 2) do θ
  sqrt(2)* [ 1-θ[1], 10*(θ[2]-θ[1]^2) ]
end

# choose the method
conf = LevenbergMarquardt_Conf()

# initial value for θ
θ_init = zeros(2)

# solve it
result = solve(nls, θ_init, conf)

# use the returned result
@assert converged(result)
θ_solution = solution(result)
```

Now to solve the same problem, but with bound constraints, proceed as
follows:

```@example session
# choose the method
conf = LevenbergMarquardt_BC_Conf()

# define bound constraints
θl = Float64[2,2]
θu = Float64[4,4]

bc = BoundConstraints(θl,θu)

# solve it
result = solve(nls, θ_init, bc, conf)

# use the returned result
@assert converged(result)
θ_solution = solution(result)
```

For furthers details see:

- [`solve(nls::AbstractNLS, θ_init::AbstractVector, conf::Abstract_Solver_Conf)`](@ref).
- [`solve(nls::AbstractNLS, θ_init::AbstractVector,bc::BoundConstraints, conf::Abstract_Solver_Conf)`](@ref).



