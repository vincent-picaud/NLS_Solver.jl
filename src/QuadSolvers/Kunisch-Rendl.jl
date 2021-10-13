using LinearAlgebra: Symmetric



@enum(BoundConstraintState_Enum,
      BoundConstraintState_LB = -1,
      BoundConstraintState_INACTIVE = 0,
      BoundConstraintState_UB = +1)

"""
With symmetric Q, modify the system

Qx+q=0

Such that

Z[i]=active_lb ⇒  x[i]=a
Z[i]=active_ub ⇒  x[i]=b

The procedure perform in place modification of `Q` and `q`
"""
function restrict_to_inactive!(Q::Symmetric,
                               q::AbstractVector,
                               Z::AbstractVector{BoundConstraintState_Enum},
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
            
            if Z[k]!=BoundConstraintState_INACTIVE

                constrained_x_k = ifelse(Z[k]==BoundConstraintState_LB,
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
            
            if Z[k]!=BoundConstraintState_INACTIVE

                constrained_x_k = ifelse(Z[k]==BoundConstraintState_LB,
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
                               Z::AbstractVector{BoundConstraintState_Enum},
                               bc::BoundConstraints{<:Real,1})
    restrict_to_inactive!(Q,q,Z,lower_bound(bc),upper_bound(bc))
end


""" 

This function updates `Z` according to `x`, `τ` and bounds `lb`, `ub`
values.

It also count how many changes have be done during this update.

No change means that the algorithm has converged.

Note: this function only modifies `Z`
"""
function update_Z!(x::AbstractVector,
                   τ::AbstractVector,
                   Z::AbstractVector{BoundConstraintState_Enum},
                   lb::AbstractVector,
                   ub::AbstractVector)
    n = length(x)

    @assert n == length(τ)
    @assert n == length(Z)
    @assert n == length(lb)
    @assert n == length(ub)
    
    count_bad_hypothesis::UInt = 0
    
    for i in 1:n

        if Z[i]==BoundConstraintState_INACTIVE

            if x[i]<=lb[i]

                count_bad_hypothesis+=1
                Z[i]=BoundConstraintState_LB
                
            elseif x[i]>=ub[i]
                
                count_bad_hypothesis+=1
                Z[i]=BoundConstraintState_UB
                
            end
            
        elseif Z[i]==BoundConstraintState_LB
            
            @assert x[i]==lb[i] "Internal error $(x[i]) != $(lb[i])"
            
            if τ[i]>0
                count_bad_hypothesis+=1
                Z[i]=BoundConstraintState_INACTIVE
            end

        else
            
            @assert x[i]==ub[i] "Internal error $(x[i]) != $(ub[i])"
            
            if τ[i]<0
                count_bad_hypothesis+=1
                Z[i]=BoundConstraintState_INACTIVE
            end
            
        end 
        
    end #  for i in 1:n

    return count_bad_hypothesis
end
function update_Z!(x::AbstractVector,
                   τ::AbstractVector,
                   Z::AbstractVector{BoundConstraintState_Enum},
                   bc::BoundConstraints{<:Real,1})
    update_Z!(x,τ,Z,lower_bound(bc),upper_bound(bc))
end 

"""
Create (x,Z) from initial guess `x_init` and bound constraints `bc`
"""
function initialize_x_Z(x_init::AbstractArray,
                        bc::BoundConstraints)

    @assert size(x_init) == size(bc)

    x = similar(x_init)
    Z = similar(Array{BoundConstraintState_Enum},axes(bc))
    lb = lower_bound(bc)
    ub = upper_bound(bc)

    for (i,(lb_i,x_init_i,ub_i)) in enumerate(zip(lb,x_init,ub))

        if x_init_i<lb_i
            x[i]=lb_i
            Z[i]=BoundConstraintState_LB
            continue
        end

        if x_init_i>ub_i
            x[i]=ub_i
            Z[i]=BoundConstraintState_UB
            continue
        end

        x[i]=x_init_i
        Z[i]=BoundConstraintState_INACTIVE
    end

    x, Z
end    

"""

Update x value according to Z
"""
function update_x!(x::AbstractArray,
                   Z::AbstractArray{BoundConstraintState_Enum},
                   bc::BoundConstraints)
    
    @assert size(x) == size(bc) == size(Z)

    for (i,(Z_i,lb_i,ub_i)) in enumerate(zip(Z,lower_bound(bc),upper_bound(bc)))
        if Z_i==BoundConstraintState_LB
            x[i]=lb_i
            continue
        end
        
        if Z_i==BoundConstraintState_UB
            x[i]=ub_i
            continue
        end
    end

    x
end 

"""
Return penalization schedule f. 

When burning_phase = false, damping factor is assumed to be one

f is a decreasing function of the iterations such that:
- f(1)=c0
- f(k)=1  if k >= k0
"""
function create_damping_schedule_exp(c0::Float64,k0::Int)
    @assert k0≥1
    @assert c0≥1

    α::Float64 = ifelse(k0>0,-log(c0)/(k0-1),zero(Float64))
    
    function damping(iter::Int)
        @assert iter≥1

        burning_phase::Bool = false
        damping_factor::Float64 = 1
        
        if iter < k0
            burning_phase = true
            damping_factor = c0*exp(α*(iter-1))
        end

        (burning_phase,damping_factor)
    end
end
function create_damping_schedule_nothing()
    
    function damping(iter::Int)
        @assert iter≥1

        (false,one(Float64))
    end 
end

"""
Put all together 
"""
function Kunisch_Rendl(Q::Symmetric{<:Real},
                       q::AbstractVector{<:Real},
                       x_init::AbstractVector{<:Real},
                       bc::BoundConstraints{<:Real,1},
                       maxIter::Int,
                       damping::Function)

    n = length(q)

    @assert (n,n) == size(Q)
    @assert n == length(bc)
    @assert first(damping(maxIter)) == false # damping stop before max iter
    
    (x, Z) = initialize_x_Z(x_init,bc)

    τ = similar(q)
    Q_tilde = similar(Q)
    q_tilde = similar(q)

    has_CV::Bool = false
    for iter in 1:maxIter
        Q_tilde .= Q
        q_tilde .= q

        (burning_phase,damping_factor) = damping(iter)

        # if still in burning phase scale diagonal 
        if burning_phase
            diagonal_indices = diagind(Q_tilde)
            Q_tilde[diagonal_indices] .*= damping_factor
        end

        # modify Q,q according to constraints
        restrict_to_inactive!(Q_tilde,q_tilde,Z,bc)

        # x = -Q^{-1}.q
        # todo; add exception
        x = -Q_tilde\q_tilde

        # update x and compute a posteriori multipliers
        update_x!(x,Z,bc)
        τ .= -Q*x-q # for perfs use mul! instead

        # update Z and count bad hypothesis
        count_bad_choice = update_Z!(x,τ,Z,bc)

        # CV check
        if (!burning_phase)&&(count_bad_choice==0)
            # clean τ
            τ_ELT = eltype(τ)
            for (i,Z_i) in enumerate(Z)
                if Z_i == BoundConstraintState_INACTIVE
                    τ[i]=zero(τ_ELT)
                end
            end
            has_CV=true
            break
        end
    end
    (has_CV,x,τ)
end


"""
Check First-Order Conditions (see
https://wiki.mcs.anl.gov/leyffer/images/0/01/07-bndCons.pdf)

If ``x^\\star=\\arg\\min f(x), x\\in[l,u]`` then:

```math
\\partial_i f(x^\\star) = \\left\\{
\\begin{array}{ll}
\\ge 0, & \\mbox{if} x^star_i = l_i \\\\
\\eq 0, & \\mbox{if} l_i \\le x^star_i \\le u_i \\\\
\\le 0, & \\mbox{if} x^star_i = u_i 
\\end{array}
\\right.
```

This is equivalent to:
```math
x^\\star = P_{[l,u]}(x^\\star-\\nabla f(x^\\star))
```

Using the previous results, this function returns:
```math
max | x^\\star - P_{[l,u]}(x^\\star-(Q.x^\\star+q))|
```
for a local optimal quantity this must be null 
"""
function check_first_order(Q::Symmetric{<:Real},
                           q::AbstractVector{<:Real},
                           xstar::AbstractVector{<:Real},
                           bc::BoundConstraints{<:Real,1})
    v = xstar-(Q*xstar+q)
    v = project!(v,bc)

    maximum(abs.(xstar .- v))
end 
