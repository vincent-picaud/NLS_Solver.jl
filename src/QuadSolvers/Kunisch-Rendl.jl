export KunischRendlConf

using LinearAlgebra: Symmetric

# Configuration
#

"""

Kunisch-Rendl solver configuration

# Regularization 

The Kunisch method is quite sensitive to ill-conditioned Q
matrix

To tackle such kind of problems we introduce a decreasing regularization
enforcing diagonal dominance

`c` is the diagonal amplification factor >= 1 applied to the
diagonal. `c` is the value at the first iteration. `c=1` at iteration
`n≥k` and the amplification has no more effect. The decreasing
schedule is empirical and defined as follows:

    cₖ = exp(-log(c₀) TODO to fix

Note: convergence tests are performed only after `k` iterations, hence
      the algorithm perform at least `k` iterations

""" 
struct KunischRendlConf <: AbstractQuadSolverConf
    iter_max::Int
    
    # Regularization
    c::Float64
    k::Int

    function KunischRendlConf(;
                              iter_max=50,
                              c=1,
                              k=1)
        @assert c≥1
        @assert iter_max>k≥1
        new(Int(iter_max),Float64(c),Int(k))
    end
end

# Output 
#
struct KunischRendlResult <: AbstractQuadSolverResult
    converged::Bool
    iter::Int
    X_solution::AbstractVector
    τ_solution::AbstractVector
end 

# ****************************************************************
# Algorithm
# ****************************************************************
#

@enum(BoundContraintState_Enum,
      BoundContraintState_LB = -1,
      BoundContraintState_INACTIVE = 0,
      BoundContraintState_UB = +1)

"""
With symmetric Q, modify the system

Qx+q=0

Such that

Z[i]=active_lb ⇒  x[i]=a
Z[i]=active_ub ⇒  x[i]=b

The procedure perform in place modification of `Q` and `q`
"""
function restrict_to_inactive!(Q::Symmetric{T},
                               q::AbstractVector{T},
                               Z::AbstractVector{BoundContraintState_Enum},
                               lb::AbstractVector{T},
                               ub::AbstractVector{T}) where {T<:Real}
    n = length(q)

    @assert (n,n)==size(Q)
    @assert n==length(Z)
    @assert n==length(lb)
    @assert n==length(ub)
    @assert Q.uplo=='L' || Q.uplo=='U'

    if Q.uplo=='L'
        
        for k = 1:n
            
            if Z[k]!=BoundContraintState_INACTIVE

                constrained_x_k = ifelse(Z[k]==BoundContraintState_LB,
                                         lb[k],
                                         ub[k])
                
                
                q[1:(k-1)] .+= constrained_x_k*@view(Q.data[k,1:(k-1)])
                Q.data[k,1:(k-1)] .= zero(T)
                
                Q_kk = max(1,abs(Q.data[k,k]))
                q[k] = -constrained_x_k*Q_kk
                Q.data[k,k] = Q_kk
                
                q[(k+1):n] .+= constrained_x_k*@view(Q.data[(k+1):n,k])
                Q.data[(k+1):n,k] .= zero(T)
                
            end # Z[k]!=inactive
        end # k = 1:n
    else
        @assert (Q.uplo=='U') ""
        
        for k = 1:n
            
            if Z[k]!=BoundContraintState_INACTIVE

                constrained_x_k = ifelse(Z[k]==BoundContraintState_LB,
                                         lb[k],
                                         ub[k])

                
                q[1:(k-1)] .+= constrained_x_k*@view(Q.data[1:(k-1),k])
                Q.data[1:(k-1),k] .= zero(T)

                Q_kk = max(1,abs(Q.data[k,k]))
                q[k] = -constrained_x_k*Q_kk
                Q.data[k,k] = Q_kk
                
                q[k+1:n] .+= constrained_x_k*@view(Q.data[k,(k+1):n])
                Q.data[k,k+1:n] .= zero(T)

            end # Z[k]!=inactive
        end # k = 1:n
    end
end

