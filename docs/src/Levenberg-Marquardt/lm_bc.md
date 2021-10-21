# Bound contrained problems

Solve an bound constrained non-linear least squares problem as defined by
[`AbstractNLS`](@ref) using the Levenberg-Marquardt method. 

To select this method define a [`Levenberg_Marquardt_BC_Conf `](@ref)
and call the [`solve`](@ref) method.


## Configuration


```@autodocs
Modules = [NLS_Solver]
Pages = ["Levenberg-Marquardt/lm_bc.jl"]
Private = false
```

## Private implementation

We discourage direct call of this **internal** function. Use the
`solve()` method instead.

```@autodocs
Modules = [NLS_Solver]
Pages = ["Levenberg-Marquardt/lm_bc.jl"]
Public = false
```
