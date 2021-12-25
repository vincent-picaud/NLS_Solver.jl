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

This package is easy to use. It follows a basic template where we have
a `solve()` function of the form:

```julia
solve(problem, algorithm_conf)::algorithm_result
```

- **problem** is problem dependant data
- **algorithm_conf** defines the chosen algorithm
- **algorithm_result** is a structure that contains the
  result.

For the moment there are only two implemented methods. The classical
Levenberg-Marquardt method for unconstrained problems. To use this
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



