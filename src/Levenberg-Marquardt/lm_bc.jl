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
function quadratic_subproblem(H::Symmetric{<:Real},
                              ∇f::AbstractVector{<:Real},
                              θ_init::AbstractVector{<:Real},
                              bc::BoundConstraints{<:Real,1},
                              conf::AbstractQuadSolverConf,
                              damping::AbstractDynamicDampingFactor,
                              max_attempt::Int)
    @assert max_attempt ≥ 1

    # translate bound onces for all: [l-θ,u-θ]
    bc_translated = bc - θ_init
    
    # buffer to avoid several allocs
    H_μI = copy(H)
    update_H_μI = begin
        view_diag_H_μI = @view H_μI[diagind(H_μI)]
        view_diag_H = @view H[diagind(H)]
        () -> begin
            μ = get_damping_factor(damping)
            view_diag_H_μI .= view_diag_H .+ μ
        end
    end
    
    # We perform several attempts with increasing μ (this is useful
    # when using some optimizers likes Kunisch-Rendl which are not
    # necessary convergent for ill-conditioned system).

    # as they are used outside the for loop
    local quad_cv_ok = false
    local result
    for attempt in 1:max_attempt
        update_H_μI()                
        result = solve(H_μI,∇f,θ_init,bc_translated,conf)

        quad_cv_ok = converged(result)
        if quad_cv_ok 
            break
        end

        # only bindings...
        θ_init = solution(result)
        damping = update_damping_factor(damping, -1.0) 
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

@doc raw"""
```julia
Levenberg_Marquardt_BC(nls::AbstractNLS,
                       θ_init::AbstractVector,
                       bc::BoundConstraints;
                       # parameters
                       max_iter::Int=50,
                       ε_grad_inf_norm::Float64=1e-8,
                       ε_step_2_norm::Float64=1e-8,
                       # initial regularization μ0=τ.|H|
                       τ::Float64=1.0e-3,
                       # quad specific
                       quad_conf::AbstractQuadSolverConf=Kunisch_Rendl_Conf(),
                       quad_max_attempt::Int=10)
```

Implementation of a "personal" Levenberg-Marquardt method that handles
bound constraints.

"""
function Levenberg_Marquardt_BC(nls::AbstractNLS,
                                θ_init::AbstractVector,
                                bc::BoundConstraints;
                                # parameters
                                max_iter::Int=50,
                                ε_grad_inf_norm::Float64=1e-8,
                                ε_step_2_norm::Float64=1e-8,
                                # initial regularization μ0=τ.|H|
                                τ::Float64=1.0e-3,
                                # quad specific
                                quad_conf::AbstractQuadSolverConf=Kunisch_Rendl_Conf(),
                                quad_max_attempt::Int=10)
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
    θ_T = eltype(θ_init)
    θ=copy(θ_init)
    # be sure that θ is inbound
    project!(θ,bc)

    (r,J)=eval_r_J(nls,θ)
    # fobj is not really used (only when we return result) hence we do
    # not create this extra variable, but only its gradient:
    ∇fobj=eval_nls_∇fobj(r,J)

    # Compute H=J'J
    #
    H=eval_nls_∇∇fobj(J)

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
            
        norm_2_step = norm(step,2)

        if norm_2_step ≤ ε_step_2_norm*max(ε_step_2_norm,norm_2_step)
            result = LevenbergMarquardt_BC_Result(_converged=true,
                                               _iter_count=iter,
                                               _fobj=eval_nls_fobj(r),
                                               _solution=θ)

            @debug "Small step" result = result
            
            return result
        end

        # compute ρ
        # (for that must update θ_new and r_new)
        δL = compute_δL_naive(J,
                              r,
                              step)
       @assert δL>0 "$δL"
        
        @. θ_new = θ + step
        project!(θ_new,bc) # be sure 
        r_new = eval_r(nls,θ_new)       
        δf = compute_δf(r,r_new)

        ρ=δf/δL


        # Accept new point?
        #
        # -> update position and check for CV
        #
        if ρ>0
            
            @. θ = θ_new
            (r,J) = eval_r_J(nls,θ_new) # r_new was already know, but not J
            ∇fobj = eval_nls_∇fobj(r,J)
            H = eval_nls_∇∇fobj(J)
            
            inf_norm_KKT = norm(∇fobj+τ,Inf)

            #            @debug "iter=$(_fmt(iter)), |step|=$(_fmt(norm_2_step)), |KKT|=$(_fmt(inf_norm_KKT)), μ=$(_fmt(get_damping_factor(damping)))" 
            
            if inf_norm_KKT ≤ ε_grad_inf_norm
                result = LevenbergMarquardt_BC_Result(_converged=true,
                                                          _iter_count=iter,
                                                   _fobj=eval_nls_fobj(r),
                                                   _solution=θ)
                
                @debug "KKT critical point" result = result
                
                return result
            end
            
        else
            # @debug "Reject point: ρ=$(_fmt(ρ)), μ=$(_fmt(get_damping_factor(damping)))"
        end
        
        
        # In all cases (accepted or not) update damping factor μ
        #
        damping = update_damping_factor(damping,ρ)
    end

    result = LevenbergMarquardt_BC_Result(_converged=false,
                                          _iter_count=max_iter,
                                          _fobj=eval_nls_fobj(r),
                                          _solution=θ)

    @debug "Too many iterations" result = result

    result
end


# ================================================================
# Use the "Solver Conf + Solve method" framework:
# 1. define "Levenberg_Marquardt_Conf"
# 2. overwrite the "solve()" function
# ================================================================


# ----------------------------------------------------------------
# 1. define "Levenberg_Marquardt_BC_Conf"
# ----------------------------------------------------------------
#

@doc raw"""
```julia
Levenberg_Marquardt_BC_Conf()
```

Configuration parameters of the Levenberg-Marquardt with bound constraints solver
"""
mutable struct Levenberg_Marquardt_BC_Conf <: AbstractNLSBCConf
    # The structure is mutable as we will add methods such as:
    # set_max_iter().
    #
    
    # Reuse LM conf
    _lm_conf::Levenberg_Marquardt_Conf
    
    # Specific to LM_BC
    #
    _quad_max_attempt::Int
    _quad_conf::AbstractQuadSolverConf
    
    # default values
        function Levenberg_Marquardt_BC_Conf(;
                                             lm_conf::Levenberg_Marquardt_Conf=Levenberg_Marquardt_Conf(),
                                             quad_max_attempt::Int=10,
                                             quad_conf::AbstractQuadSolverConf=Kunisch_Rendl_Conf())
            
            @assert quad_max_attempt ≥ 1
            
            new(lm_conf,
                quad_max_attempt,
                quad_conf)
        end
end


# ----------------------------------------------------------------
# 2. overwrite the "solve()" function
# ----------------------------------------------------------------
#
function solve(nls::AbstractNLS,
               θ_init::AbstractVector,
               bc::BoundConstraints,
               conf::Levenberg_Marquardt_BC_Conf)

    Levenberg_Marquardt_BC(nls,θ_init,bc,

                           max_iter=conf._lm_conf._max_iter,
                           ε_grad_inf_norm=conf._lm_conf._ε_grad_inf_norm,
                           ε_step_2_norm= conf._lm_conf._ε_step_2_norm,
                           
                           τ=conf._lm_conf._τ,

                           quad_conf=conf._quad_conf,
                           quad_max_attempt=conf._quad_max_attempt)
end
