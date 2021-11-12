export NLS_ForwardDiff

using ForwardDiff: jacobian

@doc raw"""

A specialization that uses the `ForwardDiff` package to compute the Jacobian.

By comparison with [`AbstractNLS`](@ref) you only have to define these
functions:
- [`parameter_size`](@ref) : returns ``n_θ``
- [`residue_size`](@ref) : returns ``n_S``
- [`eval_r`](@ref) : computation of ``\mathbf{r}``

"""
struct NLS_ForwardDiff <: AbstractNLS
    _eval_r_function::Function
    _residue_size::Int
    _parameter_size::Int
end

parameter_size(nls::NLS_ForwardDiff) = nls._parameter_size
residue_size(nls::NLS_ForwardDiff) = nls._residue_size

function eval_r(nls::NLS_ForwardDiff,θ::AbstractVector) 
    nls._eval_r_function(θ)
end

function eval_r_J(nls::NLS_ForwardDiff, θ::AbstractVector)

    r_evaluation = θ->eval_r(nls,θ)

    r = eval_r(nls,θ)
    J = jacobian(r_evaluation, θ)

    r,J
end

