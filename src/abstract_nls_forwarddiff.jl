export NLS_ForwardDiff, create_NLS_problem_using_ForwardDiff

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

@doc raw"""
```julia
create_NLS_problem_using_ForwardDiff(r::Function;domain_image_dim::Pair{Int,Int})
```
Create a NLS problem instance by sub-typing [`AbstractNLS`](@ref) type.

`r` is a function that maps a parameter vector θ to its residue. The
Jacobian matrix is computed using the `ForwardDiff` package.

# Usage example

```julia
nls = create_NLS_problem_using_ForwardDiff(2 => 2) do θ

end
```
"""
function create_NLS_problem_using_ForwardDiff(r::Function,domain_image_dim::Pair{Int,Int})
    NLS_ForwardDiff(r,last(domain_image_dim),first(domain_image_dim))
end                         

parameter_size(nls::NLS_ForwardDiff) = nls._parameter_size
residue_size(nls::NLS_ForwardDiff) = nls._residue_size

function eval_r(nls::NLS_ForwardDiff,θ::AbstractVector) 
    nls._eval_r_function(θ)
end

function eval_r_J(nls::NLS_ForwardDiff, θ::AbstractVector{T}) where {T}
    
    r_evaluation = (r,θ)->(r.=eval_r(nls,θ))
    
    r=Vector{T}(undef,residue_size(nls))

    J = jacobian(r_evaluation, r, θ)

    r,J
end

