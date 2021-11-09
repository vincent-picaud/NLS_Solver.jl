export create_nls_using_forwarddiff

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
    @assert length(θ) == parameter_size(nls)

    r = nls._eval_r_function(θ)
    
    @assert length(r) == residue_size(nls)

    r
end

function eval_r_J(nls::NLS_ForwardDiff,
                  θ::AbstractVector) 
    r_evaluation = θ -> eval_r(nls,θ)
    r = eval_r(nls,θ)
    J = jacobian(r_evaluation, θ)

    r,J
end

function create_nls_using_forwarddiff(
    eval_r_function::Function,
    residue_size::Int,
    parameter_size::Int) 

    NLS_ForwardDiff(eval_r_function,residue_size,parameter_size)
end 
