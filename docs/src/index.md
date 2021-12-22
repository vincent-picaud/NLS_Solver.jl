```@meta
CurrentModule = NLS_Solver
```

```@setup session
using NLS_Solver
```

# Introduction

This package goal is to solve nonlinear least squares problems. It
currently supports two kind of problems:

- classical nonlinear least squares:
```math
\min\limits_{\theta} \frac{1}{2}\|r(\theta)\|^2
```
- bound constrained nonlinear least squares
```math
\min\limits_{\theta_l\le\theta\le\theta_u} \frac{1}{2}\|r(\theta)\|^2
```

The implemented method is the Levenberg-Marquardt method. The bound
constrained quadratic inner problems are solved using the
Kunisch-Rendl method.

# Getting started

Please read this package [README.org](https://github.com/vincent-picaud/NLS_Solver.jl).

This example illustrates how to solve a nonlinear least squares. The objective function is the Rosenbrock function:
```math
(1-θ_1)^2 + 100 (θ_2 - θ_1^2)^2
```


```@example session
nls = create_NLS_problem_using_ForwardDiff(2 => 2) do θ
  sqrt(2)* [ 1-θ[1], 10*(θ[2]-θ[1]^2) ]
end


conf = LevenbergMarquardt_Conf()

θ_init = zeros(2)

result = solve(nls, θ_init, conf)

@assert converged(result)

θ_solution = solution(result)
```

This example shows how to solve the same problem but with bound
constraints.

```@example session
conf = LevenbergMarquardt_BC_Conf()

θl = Float64[2,2]
θu = Float64[4,4]

bc = BoundConstraints(θl,θu)

result = solve(nls, θ_init, bc, conf)

@assert converged(result)

θ_solution = solution(result)
```


# General considerations

This package is easy to use. It follows a basic template where we have
a `solve()` function of the form:

```julia
solve(problem, algorithm_conf)::algorithm_result
```

- **problem** is problem dependant data
- **algorithm_conf** is a structure that contains algorithm
  parameters. It is also useful to define the algorithm we want to
  use.
- **algorithm_result** is a returned structure that contains the
  result.

## Classical nonlinear least squares

The problem to solve is:

```math
\min\limits_{\theta} \frac{1}{2}\|r(\theta)\|^2
```

The `solve()` function to use is [`solve(nls::AbstractNLS,
θ_init::AbstractVector, conf::Abstract_Solver_Conf)`](@ref).

## Bound constrained nonlinear least squares

The problem to solve is:

```math
\min\limits_{\theta_l\le\theta\le\theta_u} \frac{1}{2}\|r(\theta)\|^2
```

The `solve()` function to use is [`solve(nls::AbstractNLS,
θ_init::AbstractVector, bc::BoundConstraints, conf::Abstract_Solver_Conf)`](@ref).



