export AbstractNLS
export parameter_size, residue_size
export eval_r!, eval_r_J!, eval_r, eval_r_J
export eval_nls_fobj, eval_nls_∇fobj!

using LinearAlgebra: dot, mul!

@doc raw"""
```julia
abstract type AbstractNLS end 
```

Define an abstract non-linear least squares problem (NLS). In our
context such problem is essentially a differentiable function

```math
r: \theta\in\mathbb{R}^{n_θ}\mapsto r(\theta)\in\mathbb{R}^{n_S}
``` 

The objective function is

```math
\frac{1}{2}\mathbf{r}^t(θ)\mathbf{r}(θ)
``` 
where ``r(θ)∈\mathbb{R}^{n_S}`` is the residue vector. Its dimension
is ``n_S``, the number of sample. ``θ`` is the vector of parameters to
optimize. Its dimension is ``n_θ``.

Its gradient is ``J r`` where ``J`` is the ``r`` Jacobian:
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

The **both** functions return `r`
"""
eval_r!(r::AbstractVector,nls::AbstractNLS,θ::AbstractVector) = error("To implement")

"""
see [`eval_r!`](@ref)
"""
function eval_r(nls::AbstractNLS, θ::AbstractVector)
    elt = eltype(θ)
    n_S = residue_size(nls)
    r = Vector{elt}(undef,n_S)

    eval_r!(r,nls,θ) # return r
end

@doc raw""" 
    eval_r_J!(r::AbstractVector, J::AbstractMatrix, nls::AbstractNLS,θ::AbstractVector)
    eval_r_J(nls::AbstractNLS,θ::AbstractVector)

In-place evaluation of residual vector ``r`` and its Jacobian matrix ``J``.

The **both** functions return `(r,J)`
"""
eval_r_J!(r::AbstractVector, J::AbstractMatrix,nls::AbstractNLS,θ::AbstractVector) = error("To implement")

"""
see [`eval_r_J!`](@ref)
"""
function eval_r_J(nls::AbstractNLS,θ::AbstractVector) 
    elt = eltype(θ)
    n_S = residue_size(nls)
    n_θ = parameter_size(nls)
    r = Vector{elt}(undef,n_S)
    J = Matrix{elt}(undef,n_θ,n_S)

    eval_r_J!(r,J,nls,θ) # return (r,J)
end

# ================================================================
#
# Some helpers...
#
@doc raw"""
```julia
eval_nls_fobj(r::AbstractVector) -> fobj
```

Compute 

```math
\frac{1}{2}\mathbf{r}^t(θ)\mathbf{r}(θ)
```
"""
eval_nls_fobj(r::AbstractVector) = dot(r,r)/2

@doc raw"""
```julia
eval_nls_∇fobj!(∇fobj::AbstractVector,
               r::AbstractVector, J::AbstractMatrix) -> nothing
```

In-place computation of gradient 

```math
\nabla\left(\frac{1}{2}\mathbf{r}^t(θ)\mathbf{r}(θ)\right) = J r
```
"""
function eval_nls_∇fobj!(∇fobj::AbstractVector,
                         r::AbstractVector, J::AbstractMatrix)
    mul!(∇fobj,J,r,1,0)

    nothing
end
    


