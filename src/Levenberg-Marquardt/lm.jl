# Unconstrained problem

@doc raw"""
```julia
abstract type AbstractNLS end 
```
Define an abstract non-linear least squares problem
"""
abstract type AbstractNLS end 

@doc raw"""
```julia
eval_fobj(nls::AbstractNLS,θ::AbstractVector)
```

Compute 

```math
\frac{1}{2}\mathbf{r}^t(θ)\mathbf{r}(θ)
```

where ``r(θ)∈\mathbb{R}^{n_S}`` is the residue vector. Its dimension
is ``n_S``, the number of sample. ``θ`` is the vector of parameters to
optimize.

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

and the differential ``dr(θ)[.]``, represented by the Jacobian matrix
``J`` of components:

```math
J_{i,j}=\partial_j r^i,\ i\in[1,n_S],\ j\in[1,n_θ]
```

where ``n_S``, the number of sample and ``n_θ`` the number of parameters.
"""
eval_fobj_J(nls::AbstractNLS,θ::AbstractVector) = error("To implement")


# Levenberg_Marquardt(
