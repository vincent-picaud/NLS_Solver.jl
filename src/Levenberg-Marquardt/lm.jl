# Unconstrained problem
# 
export LevenbergMarquardt_Conf

using LinearAlgebra: norm, I

@doc raw"""
```julia
LevenbergMarquardt(nls::AbstractNLS,
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
function LevenbergMarquardt(nls::AbstractNLS,
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
    θ=copy(θ_init)

    (r,J) = eval_r_J(nls,θ)
    ∇fobj = eval_nls_∇fobj(r,J)

    # Check CV: |∇fobj| ≤ ϵ ?
    #
    inf_norm_∇fobj = norm(∇fobj,Inf)
    if  inf_norm_∇fobj ≤ ε_grad_inf_norm

        result = LevenbergMarquardt_Result(_converged=true,
                                         _iter_count=0,
                                         _fobj=eval_nls_fobj(r),
                                         _solution=θ)
   
        # @info "Initial point was already a soluion" result = result
        
        return result
    end

    # Compute H=J'J
    #
    H = eval_nls_∇∇fobj(J)

    # Initial μ
    #
    # (maybe add (Abstract)Damping type into conf)
    #
    damping = LM_Damping(τ * norm(H,Inf))
    
    for iter ∈ 1:max_iter
        # regularize Hessian
        #
        H_μD = H + get_μ(damping)*I

        # Newton step = -inv(H).∇f
        #
        local step
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
            
            # @info "Small step" result = result
            
            return result
        end

        # Compute  δL = L(0)-L(step)
        #
        δL = compute_δL_unconstrained(∇fobj,get_μ(damping),step)
        @assert δL > 0
        
        # Compute new θ & residue
        #
        θ_new = θ + step
        r_new = eval_r(nls,θ_new)       
   
        # Compute δfobj = fobj(θ) - fobj(θ_new)
        #
        δfobj = compute_δf(r,r_new)
            
        # compute ρ = δf/δL
        #
        ρ = δfobj/δL

        # @info "LM: iter=$(_fmt(iter)), |step|=$(_fmt(norm_2_step)), |∇f|=$(_fmt(inf_norm_∇fobj)), μ=$(_fmt(get_μ(damping)))"

        # Accept new point?
        #
        # -> update position and check for CV
        #
        if ρ>0
            θ = θ_new
            r, J = eval_r_J(nls,θ) # r_new was already know, but not J
            ∇fobj = eval_nls_∇fobj(r,J)
            H = eval_nls_∇∇fobj(J)
            
            inf_norm_∇fobj = norm(∇fobj,Inf)
            if  inf_norm_∇fobj ≤ ε_grad_inf_norm

                result = LevenbergMarquardt_Result(_converged=true,
                                                 _iter_count=iter,
                                                 _fobj=eval_nls_fobj(r),
                                                 _solution=θ_new
                                                 ) 

                # @info "Found a critical point" result = result
                
                return result
            end
        end

        # In all cases (accepted or not) update damping factor μ
        #
        damping = update_μ(damping,ρ)
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
# 1. define "LevenbergMarquardt_Conf"
# 2. overwrite the "solve()" function
# ================================================================

# ----------------------------------------------------------------
# 1. define "LevenbergMarquardt_Conf"
# ----------------------------------------------------------------
#
@doc raw"""
```julia
LevenbergMarquardt_Conf()
```

Configuration parameters of the Levenberg-Marquardt solver
"""
mutable struct LevenbergMarquardt_Conf <: Abstract_Solver_Conf
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
    function LevenbergMarquardt_Conf(;
                     max_iter::Int=1000,
                     ε_grad_inf_norm::Float64=1e-8,
                     ε_step_2_norm::Float64=1e-8,
                     
                     τ::Float64=1.0e-3)
         
        # note: parameters values are controlled directly by the
        # LevenbergMarquardt() function
        
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
               conf::LevenbergMarquardt_Conf)

    LevenbergMarquardt(nls,θ_init,
                        
                        max_iter=conf._max_iter,
                        ε_grad_inf_norm=conf._ε_grad_inf_norm,
                        ε_step_2_norm= conf._ε_step_2_norm,
                        
                        τ=conf._τ)
end
