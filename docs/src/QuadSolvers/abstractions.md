# Problem abstraction

There is one [`solve`](@ref) method to solve positive definite
quadratic problem of the form:

```math
\min\limits_{\mathbf{a}\le \mathbf{x} \le \mathbf{b}} \frac{1}{2} \mathbf{x}^t\mathbf{Q}\mathbf{x} + \mathbf{q}^t\mathbf{x}
```

The first argument of the `solve()` function is the solver
configuration. This [`AbstractQuadSolverConf`](@ref) structure defines
the solver to use and its configuration. The method returns a
[`AbstractQuadSolverResult`](@ref) that contains the algorithm
results.

## Solver Configuration

```@autodocs
Modules = [NLS_Solver]
Pages = ["QuadSolvers/abstract_quadsolver_conf.jl"]
Private = false
```

### Interface (if any)

```@autodocs
Modules = [NLS_Solver]
Order   = [:function]
Pages = ["QuadSolvers/abstract_quadsolver_conf.jl"]
Private = false
```

## Solver Result

```@autodocs
Modules = [NLS_Solver]
Order   = [:type]
Pages = ["QuadSolvers/abstract_quadsolver_result.jl"]
Private = false
```

### Interface (if any)

```@autodocs
Modules = [NLS_Solver]
Order   = [:function]
Pages = ["QuadSolvers/abstract_quadsolver_result.jl"]
Private = false
```

## Solving the problem

```@autodocs
Modules = [NLS_Solver]
Order   = [:function]
Pages = ["QuadSolvers/solve.jl"]
Private = false
```

## Not exported (or specializations)

```@autodocs
Modules = [NLS_Solver]
Pages = ["QuadSolvers/abstract_quadsolver_conf.jl",
	"QuadSolvers/abstract_quadsolver_result.jl",
	"QuadSolvers/solve.jl"]
Public = false
```
