# Unconstrained problem
#
export Levenberg_Marquardt

using LinearAlgebra: norm, I

# TODO
# abstract type AbstractNLSSolverConf end
# struct Levenberg_Marquardt <:  AbstractNLSSolverConf
# ...
# end

function Levenberg_Marquardt(nls::AbstractNLS,
                             θ_init::AbstractVector;
                             # parameters
                             max_iter::Int=50,
                             ε_grad_inf_norm::Float64=1e-8,
                             ε_step_2_norm::Float64=1e-8,
                             # initial regularization
                             τ::Float64=1.0e-3,                         
                             ν_init::Float64=2.0,
                             verbose::Bool=true)
    v = ν_init
    
    # Compute, r,J, ∇fobj=J'r
    #
    n_S, n_θ = residue_size(nls),parameter_size(nls) 
    θ=copy(θ_init)

    (r,J)=eval_r_J(nls,θ)
    # fobj is not really used (only when we return result) hence we do
    # not create this extra variable, but only its gradient:
    ∇fobj=similar(r)
    eval_nls_∇fobj!(∇fobj,r,J)

    # Check CV: |∇fobj| ≤ ϵ ?
    #
    inf_norm_∇fobj = norm(∇fobj,Inf)
    if  inf_norm_∇fobj ≤ ε_grad_inf_norm
        if verbose
            @info "Already critical point CV = ok"
        end
        
        return LevenbergMarquardt_Result(_converged=true,
                                         _iter_count=0,
                                         _fobj=eval_nls_fobj(r),
                                         _solution=θ
                                         ) 
    end

    # Compute H=J'J
    #
    H=Symmetric(Matrix{eltype(r)}(undef,n_θ,n_θ))
    eval_nls_∇∇fobj!(H,J)

    # Initial μ
    #
    μ = τ * norm(H,Inf)

    # Some buffers
    #
    H_μD = similar(H) # H + μ.I
    step = similar(θ)
    θ_new = similar(θ)
    r_new = similar(r)
    
    for iter ∈ 1:max_iter
        # regularize Hessian
        #
        H_μD .= H + μ*I

        # Newton step = -inv(H).∇f
        #
        try
            step = -H_μD\∇fobj
        catch
            @warn "Unexpected singular system"
            
            return LevenbergMarquardt_Result(_converged=false,
                                             _iter_count=iter,
                                             _fobj=eval_nls_fobj(r),
                                             _solution=θ
                                             )
        end

        # Check if step not too small -> CV
        #
        norm_2_step = norm(step,2)

        if norm_2_step ≤ ε_step_2_norm*max(ε_step_2_norm,norm_2_step)
            if verbose
                @info "step too small... (TODO: check if μ is not too high
                    before saying that cv=true)"
            end
            
            return LevenbergMarquardt_Result(_converged=true,
                                             _iter_count=iter,
                                             _fobj=eval_nls_fobj(r),
                                             _solution=θ,
                                             ) 
            return true
        end

        # Compute δL variation from the quadratic model
        # δL = L(0)-L(step)
        #    = 1/2 dot( step , μ step - grad)
        #
        δL = dot(step,μ*step-∇fobj)/2         # TODO: optimize to avoid mem alloc
        @assert δL > 0
        
        # Compute new θ & residue
        #
        @. θ_new = θ + step
        eval_r!(r_new,nls,θ_new)        # TODO: optimize to avoid mem alloc
   
        # Compute δfobj = 1/2( r^2 - r_new^2 )
        # (using  r^2 - r_new^2 = (r-r_new)*(r+r_new) which is numerically better)
        #
        δfobj = dot(r-r_new,r+r_new)/2

        # compute ρ = δf/δL
        #
        ρ = δfobj/δL

        if verbose
            println("iter $iter, |step|=$norm_2_step |∇f|=$inf_norm_∇fobj μ=$μ θ=$θ_new")
        end

        if ρ>0
            # Accept new point
            #
            @. θ = θ_new
            eval_r_J!(r,J,nls,θ_new) # r_new was already know, but not J
            eval_nls_∇fobj!(∇fobj,r,J)
            eval_nls_∇∇fobj!(H,J)
            
            inf_norm_∇fobj = norm(∇fobj,Inf)
            if  inf_norm_∇fobj ≤ ε_grad_inf_norm
                if verbose
                    @info "Got a critical point CV = ok"
                end
                
                return LevenbergMarquardt_Result(_converged=true,
                                                 _iter_count=iter,
                                                 _fobj=eval_nls_fobj(r),
                                                 _solution=θ_new
                                                 ) 
            end

            μ = μ * max(1/3,1-(2*ρ-1)^3)
            ν = ν_init
        else
            # Do not accept point,
            # more regularization
            #
            μ = μ * ν
            ν = 2 * ν 
        end 
    end

    # end of loop... not convergence
    return LevenbergMarquardt_Result(_converged=false,
                                     _iter_count=max_iter,
                                     _fobj=eval_nls_fobj(r),
                                     _solution=θ
                                     ) 
end
