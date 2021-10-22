export AbstractNLS
export parameter_size, residue_size
export eval_r!, eval_r_J!, eval_r, eval_r_J
export eval_nls_fobj, eval_nls_∇fobj!, eval_nls_∇∇fobj!

using LinearAlgebra: dot, mul!
using LinearAlgebra.BLAS: BlasFloat, syrk!

@doc raw"""
```julia
abstract type AbstractNLS end 
```

Define an abstract non-linear least squares problem (NLS). In our
context such problem is essentially a differentiable function ``r``:

```math
r: \theta\in\mathbb{R}^{n_θ}\mapsto r(\theta)\in\mathbb{R}^{n_S}
``` 
where ``r(θ)∈\mathbb{R}^{n_S}`` is the residue vector. Its dimension
is ``n_S``, the number of sample. ``θ`` is the vector of parameters to
optimize. Its dimension is ``n_θ``.

The objective function is

```math
f(θ)=\frac{1}{2}\mathbf{r}^t(θ)\mathbf{r}(θ)
``` 

Its gradient ``\nabla f`` is ``J^t r`` where ``J`` is the ``r`` Jacobian:
```math
J_{i,j}=\partial_j r^i(θ),\ i\in[1,n_S],\ j\in[1,n_θ]
```

The matrix ``J^t J`` is used to approximate the Hessian. This
approximation is better when ``\| r \| ≈ 0``.

To implement a new model, you must implement:
- [`parameter_size`](@ref)
- [`residue_size`](@ref)
- [`eval_r!`](@ref)
- [`eval_r_J!`](@ref)
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
parameter_size(nls::AbstractNLS) = error("To implement")

@doc raw"""
    sample_size(nls::AbstractNLS) 

Return the dimension ``n_S`` of the residue vector ``r``.
"""
residue_size(nls::AbstractNLS) = error("To implement")

@doc raw""" 
    eval_r!(r::AbstractVector,nls::AbstractNLS,θ::AbstractVector)
    eval_r(nls::AbstractNLS,θ::AbstractVector)

In-place evaluation of residual vector ``r``

**Both** functions return `r`

"""
eval_r!(r::AbstractVector,nls::AbstractNLS,θ::AbstractVector) = error("To implement")

@doc raw""" 
    eval_r_J!(r::AbstractVector, J::AbstractMatrix, nls::AbstractNLS,θ::AbstractVector)
    eval_r_J(nls::AbstractNLS,θ::AbstractVector)

In-place evaluation of residual vector ``r`` and its
Jacobian``n_S\times n_\theta`` matrix ``J`` representing the ``dr``
differential.

**Both** functions return `(r,J)`
"""
eval_r_J!(r::AbstractVector, J::AbstractMatrix,nls::AbstractNLS,θ::AbstractVector) = error("To implement")

# ================================================================
# (Interface) convenience functions...
# ================================================================
# They look like interface functions, but but have a default
# implementation.
#

"""
    eval_r(nls::AbstractNLS, θ::AbstractVector) -> r

A convenience function that calls [`eval_r!`](@ref), but take in charge initial creation of ``r``.

"""
function eval_r(nls::AbstractNLS, θ::AbstractVector)
    elt = eltype(θ)
    n_S = residue_size(nls)
    r = Vector{elt}(undef,n_S)

    eval_r!(r,nls,θ) # return r
end


"""
     eval_r_J(nls::AbstractNLS,θ::AbstractVector) -> (r,J)

A convenience function that calls [`eval_r_J!`](@ref), but take in
charge initial creation of ``(r,J)``.
"""
function eval_r_J(nls::AbstractNLS,θ::AbstractVector) 
    elt = eltype(θ)
    n_S = residue_size(nls)
    n_θ = parameter_size(nls)
    r = Vector{elt}(undef,n_S)
    J = Matrix{elt}(undef,n_S,n_θ)

    eval_r_J!(r,J,nls,θ) # return (r,J)
end


# ================================================================
# Extra functions with implementations
# ================================================================
#
@doc raw"""
```julia
eval_nls_fobj(r::AbstractVector) -> f(θ)
```

Compute 

```math
f(θ)=\frac{1}{2}\mathbf{r}^t(θ)\mathbf{r}(θ)
```
"""
eval_nls_fobj(r::AbstractVector) = dot(r,r)/2

@doc raw"""
```julia
qeval_nls_∇fobj!(∇fobj::AbstractVector,
               r::AbstractVector, J::AbstractMatrix) -> ∇fobj
```

In-place computation of gradient 

```math
\nabla\left(\frac{1}{2}\mathbf{r}^t(θ)\mathbf{r}(θ)\right) = J^t r
```
"""
function eval_nls_∇fobj!(∇fobj::AbstractVector,
                         r::AbstractVector, J::AbstractMatrix)

    # note: no runtime penalty forJ' (this is a *lazy* operation)
    mul!(∇fobj,J',r,1,0) 
    
    ∇fobj
end
    
@doc raw"""
```julia
eval_nls_∇∇fobj!(∇∇fobj::AbstractVector,
                 J::AbstractMatrix) -> ∇∇fobj
```

In-place computation of (approximate) Hessian
```math
\nabla^2\left(\frac{1}{2}\mathbf{r}^t(θ)\mathbf{r}(θ)\right) \approx J^t.J
```
"""
function eval_nls_∇∇fobj!(∇∇fobj::Symmetric{T}, J::AbstractMatrix{T}) where {T<:BlasFloat}
    syrk!(∇∇fobj.uplo,'T',T(1),J,T(0),∇∇fobj.data)

    ∇∇fobj
end
    


