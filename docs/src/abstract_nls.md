# NLS problem abstraction

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_nls.jl"]
Private = false
```

## Unconstrained problems


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
Pages = ["abstract_nls_bc_solve.jl"]
Private = false
```

### Configuration

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_nls_bc_conf.jl"]
Private = false
```


### Result 

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_nls_bc_result.jl"]
Public = false
```

