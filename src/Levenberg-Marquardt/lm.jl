# Unconstrained problem
# 
export Levenberg_Marquardt_Conf

using LinearAlgebra: norm, I

@doc raw"""
```julia
Levenberg_Marquardt(nls::AbstractNLS,
                    θ_init::AbstractVector;
                    # parameters
                    max_iter::Int=50,
                    ε_grad_inf_norm::Float64=1e-8,
                    ε_step_2_norm::Float64=1e-8,
                    # initial regularization
                    τ::Float64=1.0e-3)
```

Implementation of a Levenberg-Marquardt method.

"""
function Levenberg_Marquardt(nls::AbstractNLS,
                             θ_init::AbstractVector;
                             # parameters
                             max_iter::Int=50,
                             ε_grad_inf_norm::Float64=1e-8,
                             ε_step_2_norm::Float64=1e-8,
                             # initial regularization
                             τ::Float64=1.0e-3)
    # Sanity check
    #
    @assert parameter_size(nls) == length(θ_init)

    @assert max_iter > 0
    @assert ε_grad_inf_norm ≥ 0
    @assert ε_step_2_norm ≥ 0
    @assert τ > 0


    # Compute, r,J, ∇fobj=J'r
    #
    n_S, n_θ = residue_size(nls),parameter_size(nls) 
    θ_T = eltype(θ_init)
    θ=copy(θ_init)

    (r,J)=eval_r_J(nls,θ)
    # fobj is not really used (only when we return result) hence we do
    # not create this extra variable, but only its gradient:
    ∇fobj=Vector{θ_T}(undef,n_θ)
    eval_nls_∇fobj!(∇fobj,r,J)

    # Check CV: |∇fobj| ≤ ϵ ?
    #
    inf_norm_∇fobj = norm(∇fobj,Inf)
    if  inf_norm_∇fobj ≤ ε_grad_inf_norm

        result = LevenbergMarquardt_Result(_converged=true,
                                         _iter_count=0,
                                         _fobj=eval_nls_fobj(r),
                                         _solution=θ)
   
        @debug "Initial point was already a soluion" result = result
        
        return result
    end

    # Compute H=J'J
    #
    H=Symmetric(Matrix{θ_T}(undef,n_θ,n_θ))
    eval_nls_∇∇fobj!(H,J)

    # Initial μ
    #
    # (maybe add (Abstract)Damping type into conf)
    #
    damping = DynamicDampingFactor(τ * norm(H,Inf))

    # Some buffers
    #
    H_μD = similar(H) # H + μ.I
    step = similar(θ)
    θ_new = similar(θ)
    r_new = similar(r)
    
    for iter ∈ 1:max_iter
        # regularize Hessian
        #
        H_μD .= H + get_damping_factor(damping)*I

        # Newton step = -inv(H).∇f
        #
        try
            step = -H_μD\∇fobj
        catch
            result = LevenbergMarquardt_Result(_converged=false,
                                               _iter_count=iter,
                                               _fobj=eval_nls_fobj(r),
                                               _solution=θ
                                               )

            @warn "Leaving LM (singular system)" result = result

            return result
        end

        # Check if step not too small -> CV
        #
        norm_2_step = norm(step,2)

        if norm_2_step ≤ ε_step_2_norm*max(ε_step_2_norm,norm_2_step)

            result = LevenbergMarquardt_Result(_converged=true,
                                             _iter_count=iter,
                                             _fobj=eval_nls_fobj(r),
                                             _solution=θ,
                                               )
            
            @debug "Small step" result = result
            
            return result
        end

        # Compute δL variation from the quadratic model
        # δL = L(0)-L(step)
        #    = 1/2 dot( step , μ step - grad)
        #
        δL = dot(step,get_damping_factor(damping)*step-∇fobj)/2         # TODO: optimize to avoid mem alloc
        @assert δL > 0
        
        # Compute new θ & residue
        #
        @. θ_new = θ + step
        eval_r!(r_new,nls,θ_new)       
   
        # Compute δfobj = 1/2( r^2 - r_new^2 )
        # (using  r^2 - r_new^2 = (r-r_new)*(r+r_new) which is numerically better)
        #
        δfobj = dot(r-r_new,r+r_new)/2

        # compute ρ = δf/δL
        #
        ρ = δfobj/δL

        #        @debug "LM: iter=$(_fmt(iter)), |step|=$(_fmt(norm_2_step)), |∇f|=$(_fmt(inf_norm_∇fobj)), μ=$(_fmt(get_damping_factor(damping)))"

        # Accept new point?
        #
        # -> update position and check for CV
        #
        if ρ>0
            @. θ = θ_new
            eval_r_J!(r,J,nls,θ_new) # r_new was already know, but not J
            eval_nls_∇fobj!(∇fobj,r,J)
            eval_nls_∇∇fobj!(H,J)
            
            inf_norm_∇fobj = norm(∇fobj,Inf)
            if  inf_norm_∇fobj ≤ ε_grad_inf_norm

                result = LevenbergMarquardt_Result(_converged=true,
                                                 _iter_count=iter,
                                                 _fobj=eval_nls_fobj(r),
                                                 _solution=θ_new
                                                 ) 

                @debug "Found a critical point" result = result
                
                return result
            end
        end

        # In all cases (accepted or not) update damping factor μ
        #
        damping = update_damping_factor(damping,ρ)
    end

    # end of loop... not convergence
    return LevenbergMarquardt_Result(_converged=false,
                                     _iter_count=max_iter,
                                     _fobj=eval_nls_fobj(r),
                                     _solution=θ
                                     ) 
end


# ================================================================
# Use the "Solver Conf + Solve method" framework:
# 1. define "Levenberg_Marquardt_Conf"
# 2. overwrite the "solve()" function
# ================================================================

# ----------------------------------------------------------------
# 1. define "Levenberg_Marquardt_Conf"
# ----------------------------------------------------------------
#
@doc raw"""
```julia
Levenberg_Marquardt_Conf()
```

Configuration parameters of the Levenberg-Marquardt solver
"""
mutable struct Levenberg_Marquardt_Conf <: AbstractNLSConf
    # The structure is mutable as we will add methods such as:
    # set_max_iter().
    #

    # Related to CV test
    #
    _max_iter::Int
    _ε_grad_inf_norm::Float64
    _ε_step_2_norm::Float64
    
    # initial regularization μ = τ max_ij(|H_ij|)
    _τ::Float64

    # default values
    function Levenberg_Marquardt_Conf(;
                     max_iter::Int=1000,
                     ε_grad_inf_norm::Float64=1e-8,
                     ε_step_2_norm::Float64=1e-8,
                     
                     τ::Float64=1.0e-3)
         
        # note: parameters values are controlled directly by the
        # Levenberg_Marquardt() function
        
        new(max_iter,
            ε_grad_inf_norm,
            ε_step_2_norm,
            
            τ)
    end
end
# TODO: add stuff like mmax_iter(), set_max_iter()...
#

# ----------------------------------------------------------------
# 2. overwrite the "solve()" function
# ----------------------------------------------------------------
#
function solve(nls::AbstractNLS,
               θ_init::AbstractVector,
               conf::Levenberg_Marquardt_Conf)

    Levenberg_Marquardt(nls,θ_init,
                        
                        max_iter=conf._max_iter,
                        ε_grad_inf_norm=conf._ε_grad_inf_norm,
                        ε_step_2_norm= conf._ε_step_2_norm,
                        
                        τ=conf._τ)
end
