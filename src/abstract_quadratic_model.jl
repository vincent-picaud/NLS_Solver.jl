export AbstractQuadraticModel

@doc raw"""

Abstract quadratic model of the form:

```math
f(θ+δθ) = f(θ) + ⟨ ∇f, δθ ⟩ + \frac{1}{2} ⟨ ∇^2f⋅δθ, δθ ⟩
``` 
# Motivation

This is an easy "generalization" of the Levenberg-Marquardt algorihtm
which is nothing more than a trust-region method.

Introducing this extra abstraction allows us to handle some extra
nonlinear least squares problems like this one:
```math
\min\limits_{θ ≥ 0} \frac{1}{2}\| r(θ) \|_2^2 + λ \| θ \|_1
```
which be expressed as:
```math
\min\limits_{θ ≥ 0} f(θ) + ⟨ J^t r + λ 1, δθ ⟩ + \frac{1}{2} ⟨ (J^tJ)⋅δθ, δθ ⟩
```

# Interface
- [`parameter_size`](@ref) : returns ``n_θ``
- [`eval_f`](@ref) : computation of ``f``
- [`eval_f_∇f`](@ref) : computation of ``(f,∇f)``
- [`eval_f_∇f_∇∇f`](@ref) : computation of ``(f,∇f,∇∇f)``

"""
abstract type AbstractQuadraticModel end


# - [`eval_f`](@ref) : computation of ``f``
# - [`eval_f_∇f`](@ref) : computation of ``(f,∇f)``
# - [`eval_f_∇f_∇∇f`](@ref) : computation of ``(f,∇f,∇∇f)``
# + TODO wrapper
