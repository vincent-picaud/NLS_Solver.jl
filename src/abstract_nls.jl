export AbstractNLS
export parameter_size, residue_size
export eval_r, eval_r_J
export eval_nls_fobj, eval_nls_∇fobj, eval_nls_∇∇fobj

using LinearAlgebra: dot, mul!
using LinearAlgebra.BLAS: BlasFloat, syrk!, gemv!

@doc raw"""
```julia
abstract type AbstractNLS end 
```

Defines an abstract non-linear least squares problem (NLS). In our
context such problem is essentially a differentiable function ``\mathbf{r}``:

```math
\mathbf{r}: \theta\in\mathbb{R}^{n_θ}\mapsto \mathbf{r}(\mathbf{\theta})\in\mathbb{R}^{n_S}
``` 
where:
- ``\mathbf{r}(\mathbf{θ})∈\mathbb{R}^{n_S}`` is the residue vector,
- ``\mathbf{θ}∈\mathbb{R}^{n_θ}`` is the parameter vector to be optimized

The objective function to minimize is:

```math
f(θ)=\frac{1}{2}\| \mathbf{r}(θ) \|^2
``` 

The classical approach uses a linear approximation of ``\mathbf{r}``:
```math
\mathbf{r}(\mathbf{θ}+δ\mathbf{θ})\approx \mathbf{r}(\mathbf{θ}) + \mathbf{J}(\mathbf{θ})\cdot δ\mathbf{θ}
``` 
where ``\mathbf{J}`` is the Jacobian:
```math
\mathbf{J}_{i,j}=\partial_j r^i(\mathbf{θ}),\ i\in[1,n_S],\ j\in[1,n_θ]
```
This leads to
```math
f(\mathbf{θ}+δ\mathbf{θ})\approx f(\mathbf{θ}) + \langle \nabla f, δ\mathbf{θ} \rangle + \frac{1}{2}  \langle \nabla^2 f \cdot δ\mathbf{θ},  δ\mathbf{θ} \rangle
```

Where the gradient ``\nabla f`` is ``\mathbf{J}^t \mathbf{r}`` and the
(approximate) Hessian ``\nabla^2 f`` is ``\mathbf{J}^t \mathbf{J}``.

To implement such model, you must define the following functions:
- [`parameter_size`](@ref) : returns ``n_θ``
- [`residue_size`](@ref) : returns ``n_S``
- [`eval_r`](@ref) : compute ``\mathbf{r}``
- [`eval_r_J`](@ref) : compute ``(\mathbf{r}, \mathbf{J})``
"""
abstract type AbstractNLS end 

# ================================================================
# Interface...
# ================================================================
#
@doc raw"""
    parameter_size(nls::AbstractNLS) 

Return the dimension ``n_θ`` of the parameter vector ``θ``.
"""
parameter_size(nls::AbstractNLS) = @assert(false,"To implement")

@doc raw"""
    sample_size(nls::AbstractNLS) 

Return the dimension ``n_S`` of the residue vector ``r``.
"""
residue_size(nls::AbstractNLS) = @assert(false,"To implement")

@doc raw""" 
```julia
eval_r(nls::AbstractNLS,
        θ::AbstractVector) -> r
```

Compte the residual vector ``\mathbf{r}``
"""
eval_r(nls::AbstractNLS,
       θ::AbstractVector) = @assert(false,"To implement")

@doc raw""" 
```julia
eval_r_J(nls::AbstractNLS,θ::AbstractVector) -> (r,J)
```

Compute the residual the vector ``\mathbf{r}`` and its Jacobian ``\mathbf{J}`` 
"""
eval_r_J(nls::AbstractNLS,
         θ::AbstractVector) = @assert(false,"To implement")

# ================================================================
# Convenience functions...
# ================================================================
#

@doc raw"""
```julia
eval_nls_fobj(r::AbstractVector{T}) -> f(θ)
```

Compute ``f(θ)=\frac{1}{2}\| \mathbf{r}(\mathbf{θ}) \|^2``
"""
eval_nls_fobj(r::AbstractVector) = dot(r,r)/2

@doc raw"""
```julia
eval_nls_∇fobj(r,J) -> ∇fobj
```

Compute the gradient: ``\nabla f(\mathbf{θ}) = \mathbf{J}^t\mathbf{r}``
"""
function eval_nls_∇fobj(r::AbstractVector, J::AbstractMatrix)
    J'*r
end
    
@doc raw"""
```julia
eval_nls_∇∇fobj(J) -> ∇∇fobj
```

Compute the (approximate) Hessian: ``\nabla^2 f(\mathbf{θ}) = \mathbf{J}^t\mathbf{J}``
"""
function eval_nls_∇∇fobj(J::AbstractMatrix) 
    Symmetric(J'*J)
end
    


