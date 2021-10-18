# Unconstrained problem
#
export Levenberg_Marquardt

using LinearAlgebra: norm

# TODO
# abstract type AbstractNLSSolverConf end
# struct Levenberg_Marquardt <:  AbstractNLSSolverConf
# ...
# end

function Levenberg_Marquardt(nls::AbstractNLS,
                             θ_init::AbstractVector;
                             # parameters
                             max_iter::Int=50,
                             ε_grad_inf_norm::Float64=1e-6,
                             ε_step_2_norm::Float64=1e-6,
                             # initial regularization
                             τ::Float64=1.0e-3,                         
                             ν::Float64=2.0,
                             verbose::Bool=true)
    # Compute, r,J, ∇fobj=J'r
    #
    n_S, n_θ = residue_size(nls),parameter_size(nls) 
    θ=copy(θ_init)

    (r,J)=eval_r_J(nls,θ)
    ∇fobj=similar(r)
    eval_nls_∇fobj!(∇fobj,r,J)

    # Check CV: |∇fobj| ≤ ϵ ?
    #
    inf_norm_∇fobj = norm(∇fobj,Inf)
    if  inf_norm_∇fobj ≤ ε_grad_inf_norm
        @info "Already critical point CV = ok"
        return true
    end

    # Compute H=J'J
    #
    H=Symmetric(Matrix{eltype(r)}(undef,n_θ,n_θ))
    eval_nls_∇∇fobj!(H,J)

    true
end
