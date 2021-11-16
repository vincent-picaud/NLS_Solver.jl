# API index 

The global API index is as follows:

```@index
Modules = [NLS_Solver]
```

## Misc

This section contains documented stuff that still have not been dispatched
across the other documentation pages. 

**CAVEAT:** some stuff can be documented without begin exported.

### Public

```@autodocs
Modules = [NLS_Solver]
Pages = ["QuadSolvers/boundconstraint_enum.jl",
         "QuadSolvers/misc.jl",
		 "Levenberg-Marquardt/damping.jl",
		 "Levenberg-Marquardt/lm_result.jl",
		 ]
Private = false	
```

### Private

```@autodocs
Modules = [NLS_Solver]
Pages = ["QuadSolvers/boundconstraint_enum.jl",
         "QuadSolvers/misc.jl",
		 "Levenberg-Marquardt/damping.jl",
		 "Levenberg-Marquardt/rho.jl",
		 "Levenberg-Marquardt/lm_result.jl",
		 ]
Public = false	
```
