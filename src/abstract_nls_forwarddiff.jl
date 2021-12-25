export NLS_ForwardDiff, create_NLS_problem_using_ForwardDiff

using ForwardDiff: jacobian

@doc raw"""

```julia
struct NLS_ForwardDiff <: AbstractNLS
    ...
end

```

A specialization that uses the `ForwardDiff` package to compute the Jacobian.

By comparison with [`AbstractNLS`](@ref) you only have to define these
functions:
- [`parameter_size`](@ref) : returns ``n_θ``
- [`residue_size`](@ref) : returns ``n_S``
- [`eval_r`](@ref) : computation of ``\mathbf{r}``

See: [`create_NLS_problem_using_ForwardDiff`](@ref) 
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

Creates an [`AbstractNLS`](@ref) specialized instance where the
[`eval_r_J`](@ref) function is automatically defined using automatic
differentiation.

- `r` is a function that maps a parameter vector θ to its residue. The
  Jacobian matrix is computed using the `ForwardDiff` package.
- `domain_image_dim` is a pair of the form `θ length => r length` that
  defines domain and codomain dimensions.

# Usage example

An example defining the Rosenbrock function

```math
\frac{1}{2}\|r(\theta)\|^2\text{ where }r = \sqrt{2} \left( \begin{array}{c}  1-\theta_1 \\ 10(\theta_2-\theta_1^2) \end{array} \right)
```

```julia
nls = create_NLS_problem_using_ForwardDiff(2 => 2) do θ
     sqrt(2) sqrt(2)* [ 1-θ[1], 10*(θ[2]-θ[1]^2) ]* [ 1-θ[1], 10*(θ[2]-θ[1]^2) ]
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

