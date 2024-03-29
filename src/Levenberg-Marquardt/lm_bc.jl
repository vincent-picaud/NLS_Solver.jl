# Problems with bound constraints
#
export LevenbergMarquardt_BC_Conf
export set_max_iteration!, set_ε_grad_Inf_norm!, set_ε_step_Inf_norm! 

using LinearAlgebra: norm, I

#
# in-place diagonal update ``diag(H_μI) = diag(H) + μ
#
function diagonal_update!(H_μI::AbstractMatrix,H::AbstractMatrix,μ::Number)
    view_diag_H_μI = @view H_μI[diagind(H_μI)]
    view_diag_H = @view H[diagind(H)]

    view_diag_H_μI .= view_diag_H .+ μ
end

@doc raw"""

Solve
```math
\min\limits_{h∈[θ^l-θ,θ^u-θ]}\frac{1}{2}h^t.(H+μI).h + ∇f^t.h
```

We use the quadratic model of ``f``, the bound contraints are such
that the step ``h`` makes the update ``x+h`` falls in the ``[θ^l,θ^u]`` bound.

"""
function quadratic_subproblem(H::Symmetric{<:Real},
                              ∇f::AbstractVector{<:Real},
                              θ_init::AbstractVector{<:Real},
                              bc::BoundConstraints{<:Real,1},
                              conf::Abstract_BC_QuadSolver_Conf,
                              damping::LM_Damping,
                              max_attempt::Int)
    @assert max_attempt ≥ 1

    # translate bound onces for all: [l-θ,u-θ]
    bc_translated = bc - θ_init
    
    # buffer to avoid several allocs
    #
    # note: to support StaticsArrays
    #
    #       H_μI = copy(H)
    #
    #       is not the right implementation, one must use similar
    #       first (in order to get a MMatrix from a SMatrix).
    #
    H_μI = similar(H)
    H_μI .= H
    
    # We perform several attempts with increasing μ (this is useful
    # when using some optimizers likes Kunisch-Rendl which are not
    # necessary convergent for ill-conditioned system).

    # as they are used outside the for loop
    local quad_cv_ok = false
    local result

    for attempt in 1:max_attempt
        
        diagonal_update!(H_μI,H,get_μ(damping))

        result = solve(H_μI,∇f,θ_init,bc_translated,conf)

        quad_cv_ok = converged(result)
        if quad_cv_ok 
            break
        end

        # only bindings...
        θ_init = solution(result)
        damping = increase_μ(damping)
    end

    if !quad_cv_ok
        @warn "Inner quad solver did not cv after $max_attempt attempts"
    end
    
    (result,damping)
end


# Note: Convergence check replace Euler CN |∇f| ≤ ϵ, we use the
# condition: max | x-P[a,b](x-∇f) | as in check_first_order
#
function LevenbergMarquardt_BC(nls::AbstractNLS,
                                θ_init::AbstractVector,
                                bc::BoundConstraints;
                                # parameters
                                max_iter::Int=50,
                                ε_grad_Inf_norm::Float64=1e-8,
                                ε_step_Inf_norm::Float64=1e-8,
                                # initial regularization μ0=τ.|H|
                                τ::Float64=1.0e-3,
                                # quad specific
                                quad_conf::Abstract_BC_QuadSolver_Conf=Kunisch_Rendl_Conf(),
                                quad_max_attempt::Int=10)
    # Sanity check
    #
    @assert length(bc) == length(θ_init)
    @assert parameter_size(nls) == length(θ_init)
    
    @assert max_iter > 0
    @assert ε_grad_Inf_norm ≥ 0
    @assert ε_step_Inf_norm ≥ 0
    @assert τ > 0

    # Initialization
    #
    θ=copy(θ_init)

    # be sure that θ is inbound
    project!(θ,bc)

    (r,J) = eval_r_J(nls,θ)
    ∇fobj = eval_nls_∇fobj(r,J)
    H = eval_nls_∇∇fobj(J)

    # Initial μ
    #
    # (maybe add (Abstract)Damping type into conf)
    #
    damping = LM_Damping(τ * norm(H,Inf))

    # Some buffers
    #
    for iter ∈ 1:max_iter
        # Find h by solving a bound constrained quadratic problem
        #
        (quad_result, damping) = quadratic_subproblem(H,
                                                      ∇fobj,
                                                      θ,
                                                      bc,
                                                      quad_conf,
                                                      damping,
                                                      quad_max_attempt) 

        if !converged(quad_result)
            @warn "LM_BC: cannot solve inner quadratic problem... Abort..."

            # we reuse last valid data (θ, ∇fobj...) and do not try to
            # use those in the quad_result structure.
            #
            return LevenbergMarquardt_BC_Result(_converged=false,
                                             _iter_count=iter,
                                             _fobj=eval_nls_fobj(r),
                                             _solution=θ)           
        end

        # Update step
        #
        step = solution(quad_result)
        τ = multiplier_τ(quad_result)
            
        norm_Inf_step = norm(step,Inf)

        if norm_Inf_step ≤ ε_step_Inf_norm
            result = LevenbergMarquardt_BC_Result(_converged=true,
                                                  _iter_count=iter,
                                                  _fobj=eval_nls_fobj(r),
                                                  _solution=θ)

            # @info "Small step" result = result
            
            return result
        end

        # compute ρ
        # (for that must update θ_new and r_new)
        #
        δL = compute_δL_constrained(∇fobj,
                                    get_μ(damping),
                                    τ,
                                    step)
        
        # δL = compute_δL_naive(J,
        #                       r,
        #                       step)
        @assert δL>0 "$δL"
        
        θ_new = θ + step
        project!(θ_new,bc) # Avoid small numerical errors that can
                           # violate bound constraints
        r_new = eval_r(nls,θ_new)
        δf = compute_δf(r,r_new)

        ρ=δf/δL


        # Accept new point?
        #
        # -> update position and check for CV
        #
        if ρ>0
            
            θ = θ_new
            r, J = eval_r_J(nls,θ) # r_new was already know, but not J
            ∇fobj = eval_nls_∇fobj(r,J)
            H = eval_nls_∇∇fobj(J)
            
            inf_norm_KKT = norm(∇fobj+τ,Inf)

            # @info "iter=$(_fmt(iter)), |step|=$(_fmt(norm_Inf_step)), |KKT|=$(_fmt(inf_norm_KKT)), μ=$(_fmt(get_μ(damping)))" 
            
            if inf_norm_KKT ≤ ε_grad_Inf_norm
                result = LevenbergMarquardt_BC_Result(_converged=true,
                                                      _iter_count=iter,
                                                      _fobj=eval_nls_fobj(r),
                                                      _solution=θ)
                
                # @info "KKT critical point" result = result
                
                return result
            end
            
        else
            # @info "Reject point: ρ=$(_fmt(ρ)), μ=$(_fmt(get_μ(damping)))"
        end
        
        
        # In all cases (accepted or not) update damping factor μ
        #
        damping = update_μ(damping,ρ)
    end

    result = LevenbergMarquardt_BC_Result(_converged=false,
                                          _iter_count=max_iter,
                                          _fobj=eval_nls_fobj(r),
                                          _solution=θ)

    @warn "Too many iterations" result = result

    result
end


# ================================================================
# Use the "Solver Conf + Solve method" framework:
# 1. define "LevenbergMarquardt_Conf"
# 2. overwrite the "solve()" function
# ================================================================


# ----------------------------------------------------------------
# 1. define "LevenbergMarquardt_BC_Conf"
# ----------------------------------------------------------------
#

@doc raw"""
```julia
mutable struct LevenbergMarquardt_BC_Conf <: Abstract_Solver_Conf
    ...
end
```

Use this constructor
```julia
LevenbergMarquardt_BC_Conf()
```
to initialize the bound constrained Levenberg-Marquardt solver default configuration
parameters.

To solve a problem with this method, you must then call 
[`solve(nls::AbstractNLS, θ_init::AbstractVector, bc::BoundConstraints, conf::Abstract_BC_Solver_Conf)`](@ref) 

See: 
- [`set_max_iteration!(conf::LevenbergMarquardt_BC_Conf,max_iter::Int)`](@ref) 
- [`set_ε_grad_Inf_norm!(conf::LevenbergMarquardt_BC_Conf,ε_grad_Inf_norm::Float64)`](@ref) 
- [`set_ε_step_Inf_norm!(conf::LevenbergMarquardt_BC_Conf,ε_step_Inf_norm::Float64)`](@ref) 

"""
mutable struct LevenbergMarquardt_BC_Conf <: Abstract_BC_Solver_Conf
    # The structure is mutable as we will add methods such as:
    # set_max_iter().
    #
    
    # Reuse LM conf
    _lm_conf::LevenbergMarquardt_Conf
    
    # Specific to LM_BC
    #
    _quad_max_attempt::Int
    _quad_conf::Abstract_BC_QuadSolver_Conf
    
    # default values
    function LevenbergMarquardt_BC_Conf(;
                                        lm_conf::LevenbergMarquardt_Conf=LevenbergMarquardt_Conf(),
                                        quad_max_attempt::Int=10,
                                        quad_conf::Abstract_BC_QuadSolver_Conf=Kunisch_Rendl_Conf())
        
        @assert quad_max_attempt ≥ 1
        
        new(lm_conf,
            quad_max_attempt,
            quad_conf)
    end
end


@doc raw"""
```julia
set_max_iteration!(conf::LevenbergMarquardt_BC_Conf,
                   max_iter::Int)
```

Modify the maximum number of iterations

See: [`LevenbergMarquardt_BC_Conf`](@ref) 
"""
function set_max_iteration!(conf::LevenbergMarquardt_BC_Conf,max_iter::Int)
    set_max_iteration!(conf._lm_conf,max_iter)
end

@doc raw"""
```julia
set_ε_grad_Inf_norm!(conf::LevenbergMarquardt_BC_Conf,
                     ε_grad_Inf_norm::Float64)
```

Modify the stopping criterion ``|\nabla f|_\infty\le\epsilon``

See: [`LevenbergMarquardt_BC_Conf`](@ref) 
"""
function set_ε_grad_Inf_norm!(conf::LevenbergMarquardt_BC_Conf,ε_grad_Inf_norm::Float64)
    set_ε_grad_Inf_norm!(conf._lm_conf,ε_grad_Inf_norm)
end

@doc raw"""
```julia
set_ε_step_Inf_norm!(conf::LevenbergMarquardt_BC_Conf,
                     ε_step_Inf_norm::Float64)
```

Modify the stopping criterion ``|\delta x|_\infty\le\epsilon``

See: [`LevenbergMarquardt_BC_Conf`](@ref) 
"""
function set_ε_step_Inf_norm!(conf::LevenbergMarquardt_BC_Conf,ε_step_Inf_norm::Float64)
    set_ε_step_Inf_norm!(conf._lm_conf,ε_step_Inf_norm)
end

# ----------------------------------------------------------------
# 2. overwrite the "solve()" function
# ----------------------------------------------------------------
#
function solve(nls::AbstractNLS,
               θ_init::AbstractVector,
               bc::BoundConstraints,
               conf::LevenbergMarquardt_BC_Conf)

    LevenbergMarquardt_BC(nls,θ_init,bc,

                           max_iter=conf._lm_conf._max_iter,
                           ε_grad_Inf_norm=conf._lm_conf._ε_grad_Inf_norm,
                           ε_step_Inf_norm= conf._lm_conf._ε_step_Inf_norm,
                           
                           τ=conf._lm_conf._τ,

                           quad_conf=conf._quad_conf,
                           quad_max_attempt=conf._quad_max_attempt)
end
