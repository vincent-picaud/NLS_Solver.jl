```@meta
CurrentModule = NLS_Solver
```

```@setup session
using NLS_Solver
```

# Unconstrained nonlinear least squares

In this part we see how to solve problems of the form:
```math
\min\limits_{\theta} \frac{1}{2}\|r(\theta)\|^2
```

For that we have to call the

[`solve(nls::AbstractNLS, θ_init::AbstractVector, conf::Abstract_Solver_Conf)`](@ref) 

function.

Where:
- [Problem definition](@ref rosenbrock_nls) section shows how to
  create an [`AbstractNLS`](@ref) instance,
- [Choose a solver](@ref lm_solver_cong) section shows how to create
   an [`Abstract_Solver_Conf`](@ref) instance.

You can reproduce the results below using the
`sandbox/example_Rosenbrock.jl` file.


## [Problem definition](@id rosenbrock_nls)

This example uses the classical Rosenbrock function:

```math
(\theta_1,\theta_2) \mapsto (1-\theta_1)^2 + 100(\theta_2-\theta_1^2)^2
```

which can be viewed as a nonlinear least squares problem, with:

```math
\frac{1}{2}\|r(\theta)\|^2\text{ where }r = \sqrt{2} \left( \begin{array}{c}  1-\theta_1 \\ 10(\theta_2-\theta_1^2) \end{array} \right)
```

The first task is to wrap this problem. The easiest way is to provide
the residue function ``r`` and let the
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) package
computes its Jacobian. This is the role of the
[`create_NLS_problem_using_ForwardDiff`](@ref) wrapper:

```@example session 
nls = create_NLS_problem_using_ForwardDiff(2 => 2) do θ
  sqrt(2)* [ 1-θ[1], 10*(θ[2]-θ[1]^2) ]
end
nothing # hide
```

The returned `nls` is an instance of [`AbstractNLS`](@ref). 

!!! danger 
    Do not specify a type, like `Float64`
    ```julia
    nls = create_NLS_problem_using_ForwardDiff(2 => 2) do θ
        sqrt(2)* Float64[ 1-θ[1], 10*(θ[2]-θ[1]^2) ]
    end
    ```
    This would prevent the use of dual numbers to compute the Jacobian.


!!! note 
    An alternative, if you do not want to use the `do...end` syntax, would
    be:
    ```julia
    Rosenbrock(θ::AbstractVector{T}) where T = sqrt(2)* T[ 1-θ[1], 10*(θ[2]-θ[1]^2) ]

    nls = create_NLS_problem_using_ForwardDiff(Rosenbrock,2 => 2);
    ```

!!! note 
    If you do not want to use `create_NLS_problem_using_ForwardDiff` you can directly sub-type
    `AbstractNLS` by defining all the required methods. This is
    described in [direct sub-typing of AbstractNLS](@ref nls_subtyping)

## [Choose a solver](@id lm_solver_cong)

The algorithm and its parameters are defined by sub-typing
[`Abstract_Solver_Conf`](@ref). For the moment there is only one
implementation, the classical Levenberg-Marquardt method:

```@example session
conf = LevenbergMarquardt_Conf()
nothing # hide
```
You can have a look at [`LevenbergMarquardt_Conf`](@ref) for further details.

We also need a starting point for the ``\theta`` parameter vector,
here we create a zero filled vector:

```@example session
θ_init = zeros(parameter_size(nls))
```

## The `solve()` function

To solve the problem one must call the `solve()` function. For
unconstrained problems, this function has the following prototype

```julia
function solve(nls::AbstractNLS,
               θ_init::AbstractVector,
               conf::Abstract_Solver_Conf)::Abstract_Solver_Result
```

In our case this gives

```@example session
result = solve(nls, θ_init, conf)
nothing # hide
```

## Using solver result

The `solve()` function returns a [`Abstract_Solver_Result`](@ref) sub-typed
structure that contains algorithm result.

In peculiar you can check if the method has converged: 

```@example session
@assert converged(result)
```

and get the founded solution: 

```@example session
θ_solution = solution(result)
```

## [Annex: direct sub-typing of `AbstractNLS`](@id nls_subtyping)

To define a nonlinear least squares problem, you can sub-type
`AbstractNLS`. In that case, you have to define everything, including
the function that computes the residue Jacobian.

More precisely, you have to define 4 methods:
- `parameter_size`: returns the ``θ`` parameter vector length
- `residue_size`: returns the ``r`` residue vector length
- `eval_r`: computes the residue ``r`` value
- `eval_r_J`: computes the residue ``r`` value and its Jacobian
  matrix wrt to ``\theta``.

For the Rosenbrock function this gives:

```@example session
struct Rosenbrock <: NLS_Solver.AbstractNLS
end

import NLS_Solver: parameter_size, residue_size, eval_r, eval_r_J

NLS_Solver.parameter_size(::Rosenbrock) = 2
NLS_Solver.residue_size(::Rosenbrock) = 2

function NLS_Solver.eval_r(nls::Rosenbrock,θ::AbstractVector{T}) where T
    @assert length(θ)==parameter_size(nls)

    sqrt(2)* T[ 1-θ[1], 10*(θ[2]-θ[1]^2) ]
end

function NLS_Solver.eval_r_J(nls::Rosenbrock,θ::AbstractVector{T}) where T
    @assert length(θ)==parameter_size(nls)

    r = sqrt(2)* T[ 1-θ[1], 10*(θ[2]-θ[1]^2) ]
    J = sqrt(2)* T[ -1 0; -20*θ[1] 10]

    (r,J)
end

nothing # hide
```

All the remaining parts are identical to the previously seen steps:

```@example session
nls = Rosenbrock()

conf = LevenbergMarquardt_Conf()
θ_init = zeros(parameter_size(nls))
result = solve(nls, θ_init, conf)
nothing # hide
```

The solution is hopefully identical to the one we found before:
```@example session
θ_solution = solution(result)
```

