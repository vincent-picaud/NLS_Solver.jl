# Problems with bound constraints
#
export Levenberg_Marquardt_BC
#export Levenberg_Marquardt_BC_Conf

using LinearAlgebra: norm, I

@doc raw"""

Solve
```math
\min\limits_{h∈[θ^l-θ,θ^u-θ]}\frac{1}{2}h^t.(H+μI).h + ∇f^t.h
```

We use the quadratic model of ``f``, the bound contraints are such
that the step ``h`` makes the update ``x+h`` falls in the ``[θ^l,θ^u]`` bound.

"""
function quadratic_subproblem(H::Symmetric{<:Real},
                              μ::Real,
                              ∇f::AbstractVector{<:Real},
                              θ_init::AbstractVector{<:Real},
                              bc::BoundConstraints{<:Real,1},
                              conf::AbstractQuadSolverConf)

    # [l-θ,u-θ]
    bc_translated = bc - θ_init

    # H+μI
    H_plus_μI = copy(H)
    H_plus_μI[diagind(H_plus_μI)] .+= μ

    solve(H_plus_μI,∇f,θ_init,bc_translated,conf)
end  
# Add with max iter + damping factor


# Try hard to solve the previous qudaratic problem by increasing damping factor
# 
function quadratic_subproblem(H::Symmetric{<:Real},
                              ∇f::AbstractVector{<:Real},
                              θ_init::AbstractVector{<:Real},
                              bc::BoundConstraints{<:Real,1},
                              conf::AbstractQuadSolverConf,
                              damping::AbstractDampingFactor,
                              max_attempt::Int)
    @assert max_attempt ≥ 1

    local quad_cv_ok = false
    local result
    for attempt in 1:max_attempt
        result = quadratic_subproblem(H,
                                      get_damping_factor(damping),
                                      ∇f,
                                      θ_init,
                                      bc,
                                      conf)

        quad_cv_ok = converged(result)
        if quad_cv_ok 
            break
        end

        # only bindings...
        θ_init = solution(result)
        damping = increase_damping_factor(damping)
    end

    if !quad_cv_ok
        @warn "Inner quad solver did not cv after $max_attempt attempts"
    end
    
    (result,damping)
end


@doc raw"""

Compute variation of the quadratic model ``L``.

``L:\mathbb{R}^{n_θ}\rightarrow \mathbb{R}`` is defined as follows:

```math
\frac{1}{2}r^t(θ+h)r(θ+h) ≈ L(h)
```
where

```math
L(h)=\frac{1}{2}r^t(θ)r(θ) + (J^tr)^th + \frac{1}{2}h^tJ^tJh
```
Given a step ``h`` we want to compute the variation 

```math
δL = L(0)-L(h) = -(J^tr)^th - \frac{1}{2}h^tJ^tJh 
``

For the **unconstrained** case, where ``h`` is solution of ``(J^tJ+μI).h+J^tr=0`` we can compute ``δL`` by:
```math
δL = \frac{1}{2}h^t(μh-J^tr)
```

However this formula is not valid for the constrained case where: ``(J^tJ+μI).h+J^tr+τ=0``.

TODO 
```math
δL = ...
```
For the moment we do the complete computation...
"""
function compute_δL_naive(J::AbstractMatrix,
                          r::AbstractVector,
                          h::AbstractVector)
    Jtr = J'*r
    Jtr_h = dot(Jtr,h)

    Jh = J*h
    JtJh = J'Jh
    htJtJh = dot(h,JtJh)
    
    -(Jtr_h + htJtJh/2)
end
# where h is assumed to verify (J^tJ+μI).h = - ∇f
function compute_δL_unconstrained(∇f::AbstractVector,
                                  μ::Real,
                                  h::AbstractVector)
    1/2*dot(h,μ*h-∇f)
end
function compute_δL_constrained(∇f::AbstractVector,
                                μ::Real,
                                τ::AbstractVector,
                                h::AbstractVector)
    @error "To implement!"
end

@doc raw"""

Compute true variation of the real model: ``δf = \frac{1}{2}(r^t(θ)r(θ)-r^t(θ+h)r(θ+h))``

Contrary to ``δL`` things are simpler. However a trick is to use an equivalent formulation:
```math
δf = \frac{1}{2}(r^t(θ)r(θ)-r^t(θ+h)r(θ+h)) = \frac{1}{2}(r(θ)-r(θ+h))^t(r(θ)+r(θ+h))
```
that has a better numerical behavior.
"""
function compute_δf(r::AbstractVector,
                    r_new::AbstractVector)
    dot(r-r_new,r+r_new)
end

# Note: Convergence check replace Euler CN |∇f| ≤ ϵ, we use the
# condition: max | x-P[a,b](x-∇f) | as in check_first_order
#

function Levenberg_Marquardt_BC(nls::AbstractNLS,
                                θ_init::AbstractVector,
                                bc::BoundConstraints;
                                # parameters
                                quad_conf::AbstractQuadSolverConf=Kunisch_Rendl_Conf(),
                                max_iter::Int=50,
                                ε_grad_inf_norm::Float64=1e-8,
                                ε_step_2_norm::Float64=1e-8,
                                # initial regularization μ0=τ.|H|
                                τ::Float64=1.0e-3,                         
                                verbose::Bool=true)
    if verbose
        @info "Entering Levenberg_Marquardt_BC (LM_BC) $(@__FILE__):$(@__LINE__)"
        @info "θ = $θ_init, bounds: θl=$(lower_bound(bc)), θu=$(upper_bound(bc))"
    end
    
    # Sanity check
    #
    @assert length(bc) == length(θ_init)
    @assert parameter_size(nls) == length(θ_init)
    
    @assert max_iter > 0
    @assert ε_grad_inf_norm ≥ 0
    @assert ε_step_2_norm ≥ 0
    @assert τ > 0

    # Compute, r,J, ∇fobj=J'r
    #
    n_S, n_θ = residue_size(nls),parameter_size(nls) 
    θ=copy(θ_init)
    # be sure that θ is inbound
    project!(θ,bc)

    (r,J)=eval_r_J(nls,θ)
    # fobj is not really used (only when we return result) hence we do
    # not create this extra variable, but only its gradient:
    ∇fobj=similar(r)
    eval_nls_∇fobj!(∇fobj,r,J)

    # Compute H=J'J
    #
    H=Symmetric(Matrix{eltype(r)}(undef,n_θ,n_θ))
    eval_nls_∇∇fobj!(H,J)

    # Initial μ
    #
    # (maybe add (Abstract)Damping type into conf)
    #
    damping = DynamicDampingFactor(τ * norm(H,Inf))

    # Some buffers
    #
    θ_new = similar(θ)
    r_new = similar(r)

#    local quad_result

    for iter ∈ 1:max_iter
        # note: this function requires a DampingFactor struct,
        #       but damping is a DynamicDampingFactor.
        quad_damping = DampingFactor(get_damping_factor(damping))
        (quad_result, _) = quadratic_subproblem(H,
                                                ∇fobj,
                                                θ,
                                                bc,
                                                quad_conf,
                                                quad_damping,
                                                20) # max attempt

        if !converged(quad_result)
            @warn "LM_BC: cannot solve inner quadratic problem... Abort..."
            return false # return no cv
        end

        # Update step
        #
        step = solution(quad_result)
        τ = multiplier_τ(quad_result)
            
        norm_2_step = norm(step,2)

        if norm_2_step ≤ ε_step_2_norm*max(ε_step_2_norm,norm_2_step)
            if verbose
                @info "LM_BC: converged[vanishing move], |step|_2 = $norm_2_step, μ = $(get_damping_factor(damping))"
                @info "LM_BC: found solution θ=$θ"
            end
            
            # return LevenbergMarquardt_Result(_converged=true,
            #                                  _iter_count=iter,
            #                                  _fobj=eval_nls_fobj(r),
            #                                  _solution=θ,
            #                                  ) 
            return true
        end

        # compute ρ
        # (for that must update θ_new and r_new)
        δL = compute_δL_naive(J,
                              r,
                              step)
       @assert δL>0 "$δL"
        
        @. θ_new = θ + step
        project!(θ_new,bc) # be sure 
        eval_r!(r_new,nls,θ_new)
        δf = compute_δf(r,r_new)

        ρ=δf/δL


        # Accept new point?
        #
        # -> update position and check for CV
        #
        if ρ>0
            
            @. θ = θ_new
            eval_r_J!(r,J,nls,θ_new) # r_new was already know, but not J
            eval_nls_∇fobj!(∇fobj,r,J)
            eval_nls_∇∇fobj!(H,J)
            
            inf_norm_KKT = norm(∇fobj+τ,Inf)

            # Screen output
            #
            if verbose
                @info "LM_BC: iter=$iter, |step|=$norm_2_step, |KKT|=$inf_norm_KKT, μ=$(get_damping_factor(damping))"
            end
            
            if inf_norm_KKT ≤ ε_grad_inf_norm
                if verbose
                    @info "LM_BC: converged[critical point for KKT], |KKT| = $inf_norm_KKT"
                    @info "LM_BC: found solution θ=$θ"
                end
                return true;
            end

        else
            # Screen output
            #
            if verbose
                println("LM_BC: iter=$iter, Reject point: ρ=$ρ, μ=$(get_damping_factor(damping))")
            end
        end
        

        # In all cases (accepted or not) update damping factor μ
        #
        damping = update_damping_factor(damping,ρ)
    end

    return false
end
