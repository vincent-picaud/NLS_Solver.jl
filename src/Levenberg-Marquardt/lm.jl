# Unconstrained problem
# 
export LevenbergMarquardt_Conf
export set_max_iteration!, set_ε_grad_Inf_norm!, set_ε_step_Inf_norm! 

using LinearAlgebra: norm, I

#
# Levenberg-Marquardt implementation
#
function LevenbergMarquardt(nls::AbstractNLS,
                            θ_init::AbstractVector;
                            # parameters
                            max_iter::Int=50,
                            ε_grad_Inf_norm::Float64=1e-8,
                            ε_step_Inf_norm::Float64=1e-8,
                            # initial regularization
                            τ::Float64=1.0e-3)
    # Sanity check
    #
    @assert parameter_size(nls) == length(θ_init)

    @assert max_iter > 0
    @assert ε_grad_Inf_norm ≥ 0
    @assert ε_step_Inf_norm ≥ 0
    @assert τ > 0


    # Compute, r,J, ∇fobj=J'r
    #
    θ=copy(θ_init)

    (r,J) = eval_r_J(nls,θ)
    ∇fobj = eval_nls_∇fobj(r,J)

    # Check CV: |∇fobj| ≤ ϵ ?
    #
    inf_norm_∇fobj = norm(∇fobj,Inf)
    if  inf_norm_∇fobj ≤ ε_grad_Inf_norm

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
        norm_Inf_step = norm(step,Inf)

        if norm_Inf_step ≤ ε_step_Inf_norm

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

        # @info "LM: iter=$(_fmt(iter)), |step|=$(_fmt(norm_Inf_step)), |∇f|=$(_fmt(inf_norm_∇fobj)), μ=$(_fmt(get_μ(damping)))"

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
            if  inf_norm_∇fobj ≤ ε_grad_Inf_norm

                result = LevenbergMarquardt_Result(_converged=true,
                                                 _iter_count=iter,
                                                 _fobj=eval_nls_fobj(r),
                                                 _solution=θ_new) 

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
mutable struct LevenbergMarquardt_Conf <: Abstract_Solver_Conf
    ...
end
```

Use this constructor
```julia
LevenbergMarquardt_Conf()
```
to initialize the Levenberg-Marquardt solver default configuration
parameters.

To solve a problem with this method, you must then call 
[`solve(nls::AbstractNLS, θ_init::AbstractVector, conf::Abstract_Solver_Conf)`](@ref) 

See: 
- [`set_max_iteration!(conf::LevenbergMarquardt_Conf,max_iter::Int)`](@ref) 
- [`set_ε_grad_Inf_norm!(conf::LevenbergMarquardt_Conf,ε_grad_Inf_norm::Float64)`](@ref) 
- [`set_ε_step_Inf_norm!(conf::LevenbergMarquardt_Conf,ε_step_Inf_norm::Float64)`](@ref) 

"""
mutable struct LevenbergMarquardt_Conf <: Abstract_Solver_Conf
    # The structure is mutable as we will add methods such as:
    # set_max_iter().
    #

    # Related to CV test
    #
    _max_iter::Int
    _ε_grad_Inf_norm::Float64
    _ε_step_Inf_norm::Float64
    
    # initial regularization μ = τ max_ij(|H_ij|)
    _τ::Float64

    # default values
    function LevenbergMarquardt_Conf(;
                     max_iter::Int=1000,
                     ε_grad_Inf_norm::Float64=1e-8,
                     ε_step_Inf_norm::Float64=1e-8,
                     
                     τ::Float64=1.0e-3)
         
        # note: parameters values are controlled directly by the
        # LevenbergMarquardt() function
        
        new(max_iter,
            ε_grad_Inf_norm,
            ε_step_Inf_norm,
            
            τ)
    end
end

@doc raw"""
```julia
set_max_iteration!(conf::LevenbergMarquardt_Conf,
                   max_iter::Int)
```

Modify the maximum number of iterations

See: [`LevenbergMarquardt_Conf`](@ref) 
"""
function set_max_iteration!(conf::LevenbergMarquardt_Conf,max_iter::Int)
    @assert max_iter>0
    conf._max_iter=max_iter
    conf
end

@doc raw"""
```julia
set_ε_grad_Inf_norm!(conf::LevenbergMarquardt_Conf,
                     ε_grad_Inf_norm::Float64)
```

Modify the stopping criterion ``|\nabla f|_\infty\le\epsilon``

See: [`LevenbergMarquardt_Conf`](@ref) 
"""
function set_ε_grad_Inf_norm!(conf::LevenbergMarquardt_Conf,ε_grad_Inf_norm::Float64)
    @assert ε_grad_Inf_norm>0
    conf._ε_grad_Inf_norm = ε_grad_Inf_norm
    conf
end

@doc raw"""
```julia
set_ε_step_Inf_norm!(conf::LevenbergMarquardt_Conf,
                     ε_step_Inf_norm::Float64)
```

Modify the stopping criterion ``|\delta x|_\infty\le\epsilon``

See: [`LevenbergMarquardt_Conf`](@ref) 
"""
function set_ε_step_Inf_norm!(conf::LevenbergMarquardt_Conf,ε_step_Inf_norm::Float64)
    @assert ε_step_Inf_norm>0
    conf._ε_step_Inf_norm = ε_step_Inf_norm
    conf
end


# ----------------------------------------------------------------
# 2. overwrite the "solve()" function
# ----------------------------------------------------------------
#
function solve(nls::AbstractNLS,
               θ_init::AbstractVector,
               conf::LevenbergMarquardt_Conf)

    LevenbergMarquardt(nls,θ_init,
                        
                        max_iter=conf._max_iter,
                        ε_grad_Inf_norm=conf._ε_grad_Inf_norm,
                        ε_step_Inf_norm= conf._ε_step_Inf_norm,
                        
                        τ=conf._τ)
end
