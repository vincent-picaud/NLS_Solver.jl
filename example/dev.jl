# 
#using Revise
#using Kunisch

# to add 
using LinearAlgebra
using Random: seed!

@enum(BoundContraintState_Enum,
      BoundContraintState_LB = -1,
      BoundContraintState_INACTIVE = 0,
      BoundContraintState_UB = +1)


"""

# Modify system (Q,q) 
#
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



#
# Test
#
seed!(1234)

n=10
A=[Rational{Int}(1,i+j-1) for i in 1:n, j in 1:n]

Q=Symmetric(A,:U)
q=Rational{Int}[i for i in 1:n]
Z=rand(instances(BoundContraintState_Enum),n)
lb=zeros(Rational{Int},n)
ub=ones(Rational{Int},n)


(Q,q)
restrict_to_inactive!(Q, q, Z, lb, ub)
(Q,q)

Q2=Rational{Int64}[1//1 0//1 0//1 1//4 1//5 1//6 0//1 0//1 0//1 0//1; 0//1 1//1 0//1 0//1 0//1 0//1 0//1 0//1 0//1 0//1; 0//1 0//1 1//1 0//1 0//1 0//1 0//1 0//1 0//1 0//1; 1//4 0//1 0//1 1//7 1//8 1//9 0//1 0//1 0//1 0//1; 1//5 0//1 0//1 1//8 1//9 1//10 0//1 0//1 0//1 0//1; 1//6 0//1 0//1 1//9 1//10 1//11 0//1 0//1 0//1 0//1; 0//1 0//1 0//1 0//1 0//1 0//1 1//1 0//1 0//1 0//1; 0//1 0//1 0//1 0//1 0//1 0//1 0//1 1//1 0//1 0//1; 0//1 0//1 0//1 0//1 0//1 0//1 0//1 0//1 1//1 0//1; 0//1 0//1 0//1 0//1 0//1 0//1 0//1 0//1 0//1 1//1]
