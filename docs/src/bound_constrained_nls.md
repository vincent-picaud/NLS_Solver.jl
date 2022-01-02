```@meta
CurrentModule = NLS_Solver
```

```@setup session
using NLS_Solver
```

# [Bound constrained nonlinear least squares](@id bc_nls_pb)

In this part we see how to solve problems of the form:
```math
\begin{align*}
\min\limits_{\theta} & \frac{1}{2}\|r(\theta)\|^2 \\
                     & \theta_l\le\theta\le\theta_u
\end{align*}
```

Compared to the unconstrained case, there are essentially two differences: 
- the `solve()` function has an extra parameter, `bc` that store bound constraints.
- The solver category is different and `Abstract_Solver_Conf` is
  replaced by `Abstract_BC_Solver_Conf`, (where `BC` stands for bound
  constrained).

Further details can be found there:

[`solve(nls::AbstractNLS, θ_init::AbstractVector, bc::BoundConstraints, conf::Abstract_BC_Solver_Conf)`](@ref)

As before, you can reproduce the results using the
`sandbox/example_Rosenbrock.jl` file.

## Problem definition

Identical to the [unconstrained case](@ref rosenbrock_nls):

```@example session 
nls = create_NLS_problem_using_ForwardDiff(2 => 2) do θ
  sqrt(2)* [ 1-θ[1], 10*(θ[2]-θ[1]^2) ]
end
nothing # hide
```

## Choose a solver

Algorithm and its parameters are defined by sub-typing
[`Abstract_BC_Solver_Conf`](@ref). For the moment there is only one
implementation, a bound constrained version of the classical
Levenberg-Marquardt method:

```@example session
conf = LevenbergMarquardt_BC_Conf()
nothing # hide
```

We also need a starting point for the unknown parameter vector
``\theta``. Here we start at point ``(3,3)``:

```@example session
θ_init = Float64[3,3]
nothing # hide
```

## The `solve()` function

To solve the problem one must call the `solve()` function.  For bound
constrained problems, this function has the following prototype:

```julia
function solve(nls::AbstractNLS,
               θ_init::AbstractVector,
               bc::BoundConstraints,
               conf::Abstract_Solver_Conf)::Abstract_Solver_Result
```

There is an extra `bc` argument which is used to define bound
constraints. By example if constraints are ``2 \le \theta_i \le 4`` we
can proceed as follows:

```@example session
θl = Float64[2,2]
θu = Float64[4,4]

bc = BoundConstraints(θl,θu)

result = solve(nls, θ_init, bc, conf)
nothing # hide
```

## Using solver result

The `solve()` function returns a [`Abstract_BC_Solver_Result`](@ref) sub-typed
structure that contains algorithm result.

You can check if the method has converged

```@example session
@assert converged(result)
```

get the optimal ``\theta``

```@example session
θ_solution = solution(result)
```

and its associated objective function value 

```@example session
objective_value(result)
```
