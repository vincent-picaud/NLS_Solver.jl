# NLS problem abstraction

TODO: use quadsolver abstraction as doc model.

## Problem abstraction


```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_nls.jl"]
Private = false
```

## Solver 


```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_nls_solve.jl"]
Private = false
```

### Configuration

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_nls_conf.jl"]
Private = false
```


### Result 


```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_nls_result.jl"]
Private = false
```


## Solver (Bound Constrained)


```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_nls_solve_bc.jl"]
Private = false
```

### Configuration

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_nls_conf_bc.jl"]
Private = false
```


### Result 

For the moment we still use [`AbstractNLSResult`](@ref)


