# Unconstrained problems

Solve an unconstrained non-linear least squares problem as defined by
[`AbstractNLS`](@ref) using the Levenberg-Marquardt method. 

To select this method define a [`LevenbergMarquardt_Conf`](@ref) and call the [`solve`](@ref) method.


## Configuration

Use this structure to define solvers configuration parameters.

```@autodocs
Modules = [NLS_Solver]
Pages = ["Levenberg-Marquardt/lm.jl"]
Private = false
```

## Private implementation

We discourage direct call of this **internal** function. Use the
`solve()` method instead.

```@autodocs
Modules = [NLS_Solver]
Pages = ["Levenberg-Marquardt/lm.jl"]
Public = false
```
