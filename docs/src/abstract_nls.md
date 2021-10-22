# NLS problem abstraction

TODO: use quadsolver abstraction as doc model.

## Unconstrained Problems


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

The abstract result type is private, only its interface is useful and
public.

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_nls_result.jl"]
Public = false
```

#### Interface 

The public interface 

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_nls_result.jl"]
Private = false
```


## Bound Constrained Problems


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

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_nls_bc_result.jl"]
Public = false
```

