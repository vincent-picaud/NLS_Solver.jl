# NLS problem abstraction

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_nls.jl"]
Private = false
```

## Unconstrained problems


```@autodocs
Modules = [NLS_Solver]
Pages = ["solver_interface.jl"]
Private = false
```

### Configuration

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_solver_conf.jl"]
Private = false
```


### Result 

The abstract result type is private, only its interface is useful and
public.

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_solver_result.jl"]
Public = false
```

#### Interface 

The public interface 

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_solver_result.jl"]
Private = false
```


## Bound Constrained Problems


```@autodocs
Modules = [NLS_Solver]
Pages = ["bc_solver_interface.jl"]
Private = false
```

### Configuration

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_bc_solver_conf.jl"]
Private = false
```


### Result 

```@autodocs
Modules = [NLS_Solver]
Pages = ["abstract_bc_solver_result.jl"]
Public = false
```

