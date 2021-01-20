# 
using Revise
using Kunisch

# to add 
using LinearAlgebra

# Define an abstract type for solver parameters
#
abstract type Solver_Parameters
end 

struct Kunisch_Parameters <: Solver_Parameters
    max_iterations::UInt
    # Kunisch method is not guaranteed to converge adding a diagonal shifting help a lot
    # c0 : user input >= 1
    # c_n = exp(-1/k0*log(c0) (with this choice, c_k=1)
    # (k is duration of artificial shifting)
    c0::Real
    k::UInt

    function Kunisch_Parameters(max_iterations::Int,c0::Real,k::Int)
        @assert(0<= 1 <= c0 , "Bad c0=$c0 value")
        @assert(0<= k <= max_iterations)
        new(max_iterations,c0,k)
    end
    function Kunisch_Parameters()
        Kunisch_Parameters(100,10,5)
    end
    
end 


abstract type Solver_Output
end

has_converged(::Solver_Output) = @assert("not implemented!")

struct Kunisch_Output{T} <: Solver_Output
    has_converged::Bool
    solution::Vector{T} 
end

function solve(paramers::Solver_Parameters,Q::Symmetric{T,<:AbstractMatrix{T}}) where {T<:Real}
    "toto toitpotipo"
    Kunisch_Output(true,[1:4;])
end


# demo

A=[[1 2];[2 50]]
p=Kunisch_Parameters()
r=solve(p,Symmetric(A))
