# Unconstrained problem
#


function Levenberg_Marquardt(nls::AbstractNLS,
                             θ_init::AbstractVector;
                             # parameters
                             max_iter::Int=50,
                             ε_grad_inf_norm::Float64=1e-6,
                             ε_step_2_norm::Float64=1e-6,
                             # initial regularization
                             τ::Float64=1.0e-3,
                             ν::Float64=2,
                             verbose::Bool=true
                             )
    θ=copy(θ_init)
end
