```@meta
CurrentModule = NLS_Solver
```

# Abstract types

## Solver Configuration

```@autodocs
Modules = [NLS_Solver]
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
## Solving the problem

To solve an optimization problem first create a
[`AbstractQuadSolverConf`](@ref) specialized instance, and then call
the [`solve`](@ref) function.

To create a `AbstractQuadSolverConf` instance the available solvers are:
- TODO


### Interface 

```@autodocs
Modules = [NLS_Solver]
Order   = [:function]
Pages = ["QuadSolvers/abstract_quadsolver_conf.jl",
	"QuadSolvers/abstract_quadsolver_result.jl",
	"QuadSolvers/solve.jl"]
Private = false
```
