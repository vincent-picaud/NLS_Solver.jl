# Problems with bound constraints
#
export Levenberg_Marquardt_BC_Conf

using LinearAlgebra: norm, I

@doc raw"""

Solve
```math
\min\limits_{h∈[θ^l-θ,θ^u-θ]}\frac{1}{2}h^t.(H+μI).h + ∇f^t.h
```

We use the quadratic model of ``f``, the bound contraints are such
that the step ``h`` makes the update ``x+h`` falls in the ``[θ^l,θ^u]`` bound.

"""
TODO


# Note: Convergence check replace Euler CN |∇f| ≤ ϵ, we use the
# condition: max | x-P[a,b](x-∇f) | as in check_first_order
#

function Levenberg_Marquardt_BC(nls::AbstractNLS,
                             θ_init::AbstractVector;
                             bc::BoundConstraints,
                             
                             # parameters
                             max_iter::Int=50,
                             ε_grad_inf_norm::Float64=1e-8,
                             ε_step_2_norm::Float64=1e-8,
                             # initial regularization
                             τ::Float64=1.0e-3,                         
                             ν_init::Float64=2.0,
                             verbose::Bool=true)
    # Sanity check
    #
    @assert length(bc) == length(θ_init)
    @assert parameter_size(nls) == length(θ_init)
    
    @assert max_iter > 0
    @assert ε_grad_inf_norm ≥ 0
    @assert ε_step_2_norm ≥ 0
    @assert τ > 0
    @assert ν_init > 0

    # Init 
    #
    v = ν_init
end
