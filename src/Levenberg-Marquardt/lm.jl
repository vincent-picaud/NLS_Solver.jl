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
    # Compute, r,J, ∇fobj
    #
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

    
end
