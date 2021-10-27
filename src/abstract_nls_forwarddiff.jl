export create_nls_using_forwarddiff

using ForwardDiff: jacobian!

@doc raw"""

A specialization that uses the `ForwardDiff` package to compute the Jacobian.

By comparison with [`AbstractNLS`](@ref) you only have to define these
functions:
- [`parameter_size`](@ref) : returns ``n_θ``
- [`residue_size`](@ref) : returns ``n_S``
- [`eval_r!`](@ref) : in-place computation of ``\mathbf{r}``

"""
struct NLS_ForwardDiff <: AbstractNLS
    _eval_r!_function::Function
    _residue_size::Int
    _parameter_size::Int
end
parameter_size(nls::NLS_ForwardDiff) = nls._parameter_size
residue_size(nls::NLS_ForwardDiff) = nls._residue_size

function eval_r!(r::AbstractVector,nls::NLS_ForwardDiff,θ::AbstractVector) 
    @assert length(r) == residue_size(nls)
    @assert length(θ) == parameter_size(nls)

    nls._eval_r!_function(r,θ)

    r
end

function eval_r_J!(r::AbstractVector,
                   J::AbstractMatrix,
                   nls::NLS_ForwardDiff,
                   θ::AbstractVector) 
    inplace_r_evaluation = (r,θ)->eval_r!(r,nls,θ)
    jacobian!(J,inplace_r_evaluation, r, θ)

    r,J
end

function create_nls_using_forwarddiff(
    eval_r!_function::Function,
    T::DataType,
    residue_size::Int,
    parameter_size::Int) 

    NLS_ForwardDiff(eval_r!_function,residue_size,parameter_size)
end 
