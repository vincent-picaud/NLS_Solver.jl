using LinearAlgebra: Symmetric, dot, diagind

@doc raw"""
```julia
restrict_to_inactive!(Q::Symmetric,
                      q::AbstractVector,
                      Z::AbstractVector{BoundConstraint_Enum},
                      lb::AbstractVector{<:Real},
                      ub::AbstractVector{<:Real})
```

```julia
function restrict_to_inactive!(Q::Symmetric,
                               q::AbstractVector,
                               Z::AbstractVector{BoundConstraint_Enum},
                               bc::BoundConstraints{<:Real,1})
```

In-place modification of ``(Q,q)`` that produces
``(\tilde{Q},\tilde{q})`` such that the initial optimization problem:

```math
\tilde{x} = \arg\min \frac{1}{2} x^t Q x + q^tx
```
under these constraints:
```math 
x[i] = \left\{\begin{array}{ll} 
l[i], & \text{if } Z[i] = -1 \\
u[i], & \text{if } Z[i] = +1 
\end{array}\right.
```
is transformed into this linear system:
```math
\tilde{Q}\tilde{x}+\tilde{q}=0
```
"""
function restrict_to_inactive!(Q::Symmetric,
                               q::AbstractVector,
                               Z::AbstractVector{BoundConstraint_Enum},
                               lb::AbstractVector{<:Real},
                               ub::AbstractVector{<:Real})
    Q_ET = eltype(Q)
    n = length(q)

    @assert (n,n)==size(Q)
    @assert n==length(Z)
    @assert n==length(lb)
    @assert n==length(ub)
    @assert Q.uplo=='L' || Q.uplo=='U'

    if Q.uplo=='L'
        
        for k = 1:n
            
            if Z[k]!=INACTIVE_BC

                constrained_x_k = ifelse(Z[k]==ACTIVE_LB,
                                         lb[k],
                                         ub[k])
                
                
                q[1:(k-1)] .+= constrained_x_k*@view(Q.data[k,1:(k-1)])
                Q.data[k,1:(k-1)] .= zero(Q_ET)
                
                Q_kk = max(1,abs(Q.data[k,k]))
                q[k] = -constrained_x_k*Q_kk
                Q.data[k,k] = Q_kk
                
                q[(k+1):n] .+= constrained_x_k*@view(Q.data[(k+1):n,k])
                Q.data[(k+1):n,k] .= zero(Q_ET)
                
            end # Z[k]!=inactive
        end # k = 1:n
    else
        @assert (Q.uplo=='U') ""
        
        for k = 1:n
            
            if Z[k]!=INACTIVE_BC

                constrained_x_k = ifelse(Z[k]==ACTIVE_LB,
                                         lb[k],
                                         ub[k])

                
                q[1:(k-1)] .+= constrained_x_k*@view(Q.data[1:(k-1),k])
                Q.data[1:(k-1),k] .= zero(Q_ET)

                Q_kk = max(1,abs(Q.data[k,k]))
                q[k] = -constrained_x_k*Q_kk
                Q.data[k,k] = Q_kk
                
                q[k+1:n] .+= constrained_x_k*@view(Q.data[k,(k+1):n])
                Q.data[k,k+1:n] .= zero(Q_ET)

            end # Z[k]!=inactive
        end # k = 1:n
    end
end
function restrict_to_inactive!(Q::Symmetric,
                               q::AbstractVector,
                               Z::AbstractVector{BoundConstraint_Enum},
                               bc::BoundConstraints{<:Real,1})
    restrict_to_inactive!(Q,q,Z,lower_bound(bc),upper_bound(bc))
end


@doc raw""" 
```julia
update_Z!(x::AbstractVector,
          τ::AbstractVector,
          Z::AbstractVector{BoundConstraint_Enum},
          lb::AbstractVector,
          ub::AbstractVector)
```
```julia
update_Z!(x::AbstractVector,
          τ::AbstractVector,
          Z::AbstractVector{BoundConstraint_Enum},
          bc::BoundConstraints{<:Real,1})
```
This function updates `Z` according to `x`, `τ` and bounds `lb`, `ub`
values.

It also count how many changes have be done during this update.

No change means that the algorithm has converged.

**Note:** this function only modifies `Z` and return the number of bad
hypothesis.
"""
function update_Z!(x::AbstractVector,
                   τ::AbstractVector,
                   Z::AbstractVector{BoundConstraint_Enum},
                   lb::AbstractVector,
                   ub::AbstractVector)
    n = length(x)

    @assert n == length(τ)
    @assert n == length(Z)
    @assert n == length(lb)
    @assert n == length(ub)
    
    count_bad_hypothesis::UInt = 0
    
    for i in 1:n

        if Z[i]==INACTIVE_BC

            if x[i]<=lb[i]

                count_bad_hypothesis+=1
                Z[i]=ACTIVE_LB
                
            elseif x[i]>=ub[i]
                
                count_bad_hypothesis+=1
                Z[i]=ACTIVE_UB
                
            end
            
        elseif Z[i]==ACTIVE_LB
            
            @assert x[i]==lb[i] "Internal error $(x[i]) != $(lb[i])"
            
            if τ[i]>0
                count_bad_hypothesis+=1
                Z[i]=INACTIVE_BC
            end

        else
            
            @assert x[i]==ub[i] "Internal error $(x[i]) != $(ub[i])"
            
            if τ[i]<0
                count_bad_hypothesis+=1
                Z[i]=INACTIVE_BC
            end
            
        end 
        
    end #  for i in 1:n

    return count_bad_hypothesis
end
function update_Z!(x::AbstractVector,
                   τ::AbstractVector,
                   Z::AbstractVector{BoundConstraint_Enum},
                   bc::BoundConstraints{<:Real,1})
    update_Z!(x,τ,Z,lower_bound(bc),upper_bound(bc))
end 

@doc raw"""
```julia
initialize_x_Z(x_init::AbstractArray,
               bc::BoundConstraints)
```

Create (x,Z) from initial guess `x_init` and bound constraints `bc`

`Z` is created by recording how `x_init` fulfils the bound constraints `bc`.

`x` is the projection of `x_init` on the bounded domain [l,b].
"""
function initialize_x_Z(x_init::AbstractArray,
                        bc::BoundConstraints)

    @assert size(x_init) == size(bc)

    x = similar(x_init)
    Z = similar(Array{BoundConstraint_Enum},axes(bc))
    lb = lower_bound(bc)
    ub = upper_bound(bc)

    for (i,(lb_i,x_init_i,ub_i)) in enumerate(zip(lb,x_init,ub))

        if x_init_i<lb_i
            x[i]=lb_i
            Z[i]=ACTIVE_LB
            continue
        end

        if x_init_i>ub_i
            x[i]=ub_i
            Z[i]=ACTIVE_UB
            continue
        end

        x[i]=x_init_i
        Z[i]=INACTIVE_BC
    end

    x, Z
end    

@doc raw"""
```julia
update_x!(x::AbstractArray,
          Z::AbstractArray{BoundConstraint_Enum},
          bc::BoundConstraints)
```

Update x value such that:
```math 
x[i] = \left\{\begin{array}{ll} 
l[i], & \text{if } Z[i] = -1 \\
u[i], & \text{if } Z[i] = +1 
\end{array}\right.
```
When ``Z[i]=0`` the ``x[i]`` value is unaffected.
"""
function update_x!(x::AbstractArray,
                   Z::AbstractArray{BoundConstraint_Enum},
                   bc::BoundConstraints)
    
    @assert size(x) == size(bc) == size(Z)

    for (i,(Z_i,lb_i,ub_i)) in enumerate(zip(Z,lower_bound(bc),upper_bound(bc)))
        if Z_i==ACTIVE_LB
            x[i]=lb_i
            continue
        end
        
        if Z_i==ACTIVE_UB
            x[i]=ub_i
            continue
        end
    end

    x
end 

@doc raw"""
```julia
clean_τ!(τ::AbstractArray{<:Real},              
         Z::AbstractArray{BoundConstraint_Enum})
```

By definition τ=-Qx-q. If the algorithm converged, then one must have
τ[i]=0 when the constraint is inactive.

This function updates τ by overwriting τ[i]=0 when Z[i]=inactive. 
"""
function clean_τ!(τ::AbstractArray{<:Real},
                  Z::AbstractArray{BoundConstraint_Enum})

    @assert size(τ)==size(Z)
    
    τ_ELT = eltype(τ)

    for (i,Z_i) in enumerate(Z)
        if Z_i == INACTIVE_BC
            τ[i]=zero(τ_ELT)
        end
    end
    τ
end

# ****************************************************************

# Config stucture
struct Kunisch_Rendl_Conf <: Abstract_BC_QuadSolver_Conf
    _max_iter::Int
    
    function Kunisch_Rendl_Conf(;
                                max_iter::Int=50)

        @assert max_iter>0
        
        new(max_iter)
    end
end

max_iter(conf::Kunisch_Rendl_Conf) = conf._max_iter

# ****************************************************************

# Result structure
struct Kunisch_Rendl_Result <: Abstract_BC_QuadSolver_Result
    _cv::Bool
    _iter_count::Int
    _fobj::Real
    _x::AbstractVector{<:Real}
    _τ::AbstractVector{<:Real}

    # For a reason I ignore, using Base.@kwdef leads to an error with
    # Julia 1 By consequence I do not use this macro and write the
    # constructor manually
    function  Kunisch_Rendl_Result(;
                                   _cv,
                                   _iter_count,
                                   _fobj,
                                   _x,
                                   _τ)
        new(_cv,_iter_count,_fobj,_x,_τ)
    end 

end 
converged(r::Kunisch_Rendl_Result) = r._cv
iteration_count(r::Kunisch_Rendl_Result) = r._iter_count
objective_value(r::Kunisch_Rendl_Result) = r._fobj
multiplier_τ(r::Kunisch_Rendl_Result) = ReadOnlyArray(r._τ)
solution(r::Kunisch_Rendl_Result) = ReadOnlyArray(r._x)

# ****************************************************************

"""
Put all together 
"""
function Kunisch_Rendl(Q::Symmetric{<:Real},
                       q::AbstractVector{<:Real},
                       x_init::AbstractVector{<:Real},
                       bc::BoundConstraints{<:Real,1},
                       conf::Kunisch_Rendl_Conf)
    n = length(q)

    @assert (n,n) == size(Q) "(n,n) == size(Q) $((n,n)) == $(size(Q))"
    @assert n == length(bc)

    
    (x, Z) = initialize_x_Z(x_init,bc)

    τ = similar(q)
    Q_tilde = similar(Q)
    q_tilde = similar(q)

    has_CV::Bool = false
    local iter_count
    for iter in 1:max_iter(conf)
        iter_count = iter
        Q_tilde .= Q
        q_tilde .= q

        # modify Q,q according to constraints
        restrict_to_inactive!(Q_tilde,q_tilde,Z,bc)

        # x = -Q^{-1}.q
        # TODO: add exception is Q_tilde is too badly conditioned
        try
            x = -Q_tilde\q_tilde
        catch
            @warn "Singular system... Abort..."
            break 
        end 

        # update x and compute a posteriori multipliers
        update_x!(x,Z,bc)
        τ .= -Q*x-q # for perfs use mul! instead

        # update Z and count bad hypothesis
        count_bad_choice = update_Z!(x,τ,Z,bc)

        # @info "KR: iter=$(_fmt(iter)), wrong_hypothesis=$(_fmt(count_bad_choice))"
        
        # CV check
        if count_bad_choice==0
            has_CV=true
            break
        end
    end
    # prepare quantities export

    # as τ=-Q*x-q, fobj is 1/2( -τ.x + q.x)
    #
    fobj = (-dot(τ,x) + dot(q,x))/2
    # "true" 0 for inactive constraints
    multiplier_τ = clean_τ!(τ,Z)

    result = Kunisch_Rendl_Result(
        _cv=has_CV,
        _iter_count=iter_count,
        _fobj=fobj,
        _x=x, # no risk as x is local
        _τ=τ # no risk as τ is local
    )

    # @info "Leaving Kunisch-Rendl" result = result

    result
end

# Specialize the solve method <- this method is exported
#
solve(Q::Symmetric{<:Real},
      q::AbstractVector{<:Real},
      x_init::AbstractVector{<:Real},
      bc::BoundConstraints{<:Real,1},
      conf::Kunisch_Rendl_Conf) = Kunisch_Rendl(Q,q,x_init,bc,conf)
