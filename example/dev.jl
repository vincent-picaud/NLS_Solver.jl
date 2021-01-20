# 
using Revise
using Kunisch

# to add 
using LinearAlgebra

A = [1 2; 2 50]

@edit cholesky(A)

# Define an abstract type for solver parameters
#
abstract type Solver_Parameters
end 

struct Kunisch_Parameters <: Solver_Parameters
    max_iterations::UInt
    # Kunisch method is not guaranteed to converge adding a diagonal shifting help a lot
    # c_0 : user input >= 1
    # c_n = exp(-1/k0*log(c0) (with this choice, c_k=1)
    # (k is duration of artificial shifting)
    c_0::Real
    k::UInt

    function Kunisch_Parameters(max_iterations::Int,c_0::Real,k::Int)
        @assert(0<= 1 <= c_0 , "Bad c_0=$c_0 value")
        @assert(0<= k<=max_iterations)
        new(max_iterations,c_0,k)
    end
end 
    
function solve(paramers::Solver_Parameters,Q::Symmetric{T,<:AbstractMatrix{T}}) where {T<:Real}
    "toto toitpotipo"
end
