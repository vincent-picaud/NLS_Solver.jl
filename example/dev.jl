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


""" 

After having modified ``Q, q`` using [``](@ref) one solves
``x=-Q^{-1}q``. We then compute ``\tau=-(Qx+q)`` to get a posteriori
mulitpliers.

This function update `Z` according to this inputs. It returns how many
`Z` components had to be modified. No modification means that the
multipliers are compatible with the `Z` array, the algorithm has
converged.
"""
function update_Z!(x::AbstractVector,
                   τ::AbstractVector,
                   Z::AbstractVector{BoundContraintState_Enum},
                   lb::AbstractVector,
                   ub::AbstractVector)
    n = length(x)

    @assert n == length(τ)
    @assert n == length(Z)
    @assert n == length(lb)
    @assert n == length(ub)
    
    count_bad_hypothesis::UInt = 0
    
    for i in 1:n

        if Z[i]==BoundContraintState_INACTIVE

            if x[i]<=lb[i]

                count_bad_hypothesis+=1
                Z[i]=BoundContraintState_LB
                
            elseif x[i]>=ub[i]
                
                count_bad_hypothesis+=1
                Z[i]=BoundContraintState_UB
                
            end
            
        elseif Z[i]==BoundContraintState_LB
            
            @assert  x[i]==lb[i] "Internal error $(x[i]) != $(lb[i])"
            
            if τ[i]>0
                count_bad_hypothesis+=1
                Z[i]=BoundContraintState_INACTIVE
            end

        else
            
          @assert  x[i]==ub[i] "Internal error $(x[i]) != $(ub[i])"
            
            if τ[i]<0
                count_bad_hypothesis+=1
                Z[i]=BoundContraintState_INACTIVE
            end
            
        end 
        
    end #  for i in 1:n

    return count_bad_hypothesis
end

"""
Put all together 
"""
function Kunisch_Rendl(Q::Symmetric{T},
                       q::AbtractVector{T},
                       lb::AbtractVector{T},
                       ub::AbtractVector{T},
                       maxIter::UInt = 50,
                       k0::UInt = 6,
                       c0::Float64=0.1) where {T<:Real}
    n = length(q)

    
    @assert ((n==length(lb))&&
             (n==length(ub))&&
             ((n,n)==size(Q))) "Bad dimension"

    @assert count(x->!x,0.<=ub-lb)==0 "Incoherent bounds"

    x=copy(lb)
    x=min(x,ub)
    Z=map(x->ifelse(isfinite(x),active_lb,inactive),x)

    for iter in 1:maxIter

        # [Solution]
        Q_tilde=copy(Q)
        q_tilde=copy(q)

        if iter<=k0
            mu::T=c0*norm(Q,Inf)/2^iter
            for i in 1:n Q_tilde.data[i,i]+=mu end
        end

        restrict2Active!(Q_tilde,q_tilde,Z,lb,ub)

        x=-1*(Q_tilde\q_tilde);
        tau=-1*(Q*x+q)
        # [Solution]


        count_bad_hypothesis = updateZ!(x,tau,Z,lb,ub);
        
        print("\niter $iter count $count_bad_hypothesis")
        
        if (iter>k0)&&(count_bad_hypothesis == 0)
            print("\nConverged!")
            return x,tau
        end
        
    end # iter
    
    error("\nDid not converged! \n$x")
    
end
