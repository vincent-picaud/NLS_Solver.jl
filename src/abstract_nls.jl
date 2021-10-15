export AbstractNLS
export parameter_size, eval_fobj, eval_fobj_J

@doc raw"""
```julia
abstract type AbstractNLS end 
```
Define an abstract non-linear least squares problem:

```math
\frac{1}{2}\mathbf{r}^t(θ)\mathbf{r}(θ)
``` 
where ``r(θ)∈\mathbb{R}^{n_S}`` is the residue vector. Its dimension
is ``n_S``, the number of sample. ``θ`` is the vector of parameters to
optimize. Its dimension is ``n_θ``.
"""
abstract type AbstractNLS end 


"""
    parameter_size(nls::AbstractNLS) 

Return the dimension ``n_θ`` of the parameter vector ``θ``.
"""
parameter_size(nls::AbstractNLS) = error("To implement")


@doc raw"""
```julia
eval_fobj(nls::AbstractNLS,θ::AbstractVector)
```

Compute 

```math
\frac{1}{2}\mathbf{r}^t(θ)\mathbf{r}(θ)
```
"""
eval_fobj(nls::AbstractNLS,θ::AbstractVector) = error("To implement")

@doc raw"""
```julia
eval_fobj_J(nls::AbstractNLS,θ::AbstractVector)
```

Compute
```math
\frac{1}{2}\mathbf{r}^t(θ)\mathbf{r}(θ)
```

and the differential ``dr``, represented by the Jacobian matrix
``J`` of components:

```math
J_{i,j}=\partial_j r^i(θ),\ i\in[1,n_S],\ j\in[1,n_θ]
```

where ``n_S``, the number of sample and ``n_θ`` the number of parameters.
"""
eval_fobj_J(nls::AbstractNLS,θ::AbstractVector) = error("To implement")


