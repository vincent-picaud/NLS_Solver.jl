export AbstractQuadraticModel
export parameter_size, eval_f, eval_f_∇f, eval_f_∇f_∇∇f
export NLSProblemAsQuadraticModel


using LinearAlgebra: dot, Symmetric

@doc raw"""

Abstract quadratic model of the form:

```math
f(θ+δθ) = f(θ) + ⟨ ∇f, δθ ⟩ + \frac{1}{2} ⟨ ∇^2f⋅δθ, δθ ⟩
``` 
# Motivation

This is an easy "generalization" of the Levenberg-Marquardt algorithm
which is nothing more than a trust-region method.

Introducing this extra abstraction allows us to handle some extra
nonlinear least squares problems like:
```math
\min\limits_{θ ≥ 0} \frac{1}{2}\| r(θ) \|_2^2 + λ \| θ \|_1
```
which be expressed as:
```math
\min\limits_{θ ≥ 0} f(θ) + ⟨ J^t r + λ 1, δθ ⟩ + \frac{1}{2} ⟨ (J^tJ)⋅δθ, δθ ⟩
```

# Interface
- [`parameter_size`](@ref parameter_size(qm::AbstractQuadraticModel)) : returns ``n_θ``
- [`eval_f`](@ref) : computation of ``f``
- [`eval_f_∇f`](@ref) : computation of ``(f,∇f)``
- [`eval_f_∇f_∇∇f`](@ref) : computation of ``(f,∇f,∇∇f)``

# Concrete implemenentation

- [`NLSProblemAsQuadraticModel`](@ref) : encapsulate a [`AbstractNLS`](@ref) problem

"""
abstract type AbstractQuadraticModel end

@doc raw"""
```julia
parameter_size(qm::AbstractQuadraticModel) -> Int
```

Return problem dimension, a vector length
"""
parameter_size(qm::AbstractQuadraticModel) = @assert(false,"To implement")

@doc raw"""
```julia
eval_f(qm::AbstractQuadraticModel,θ::AbstractVector) -> Real
```

Compute objective function value (a real number)
"""
eval_f(qm::AbstractQuadraticModel,θ::AbstractVector) = @assert(false,"To implement")

@doc raw"""
```julia
eval_f_∇f(qm::AbstractQuadraticModel,θ::AbstractVector) -> (Real,Vector)
```

Compute objective function value (a real number) and its gradient vector.
"""
eval_f_∇f(qm::AbstractQuadraticModel,θ::AbstractVector) = @assert(false,"To implement")

@doc raw"""
```julia
eval_f_∇f_∇∇f(qm::AbstractQuadraticModel,θ::AbstractVector) -> (Real,Vector,Matrix)
```

Compute objective function value (a real number), its gradient vector
and its Hessian symmetric matrix.

"""
eval_f_∇f_∇∇f(qm::AbstractQuadraticModel,θ::AbstractVector) = @assert(false,"To implement")

# ================================================================

@doc raw"""
```julia
NLSProblemAsQuadraticModel(nls::AbstractNLS)
```

Create a [`AbstractQuadraticModel`](@ref) from an [`AbstractNLS`](@ref)
instance.

"""
struct NLSProblemAsQuadraticModel <: AbstractQuadraticModel
    _nls::AbstractNLS
end

parameter_size(qm::NLSProblemAsQuadraticModel) = parameter_size(qm._nls)

function eval_f(qm::NLSProblemAsQuadraticModel,θ::AbstractVector)
    r = eval_r(qm._nls,θ)

    dot(r,r)/2
end

function eval_f_∇f(qm::NLSProblemAsQuadraticModel,θ::AbstractVector)
    (r,J)=eval_r_J(qm._nls,θ)

    (dot(r,r)/2,J'*r)
end

function eval_f_∇f_∇∇f(qm::NLSProblemAsQuadraticModel,θ::AbstractVector)
    (r,J)=eval_r_J(qm._nls,θ)

    (dot(r,r)/2,J'*r,Symmetric(J'*J))
end

