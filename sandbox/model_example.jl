#
# A model example
#
abstract type Abstract_Fit_Model end

parameter_size(::Abstract_Fit_Model) = @assert(false,"To implement!")

# Here eval Y_i for point X_i. When performing fit we have a collection of such X_i
#
eval_y(::Abstract_Fit_Model,X_i::Any,θ::AbstractVector)  = @assert(false,"To implement!")

# ----------------------------------------------------------------

abstract type Abstract_Fit_Model_Peak <: Abstract_Fit_Model end

# ----------------------------------------------------------------

struct Gaussian_Peak <: Abstract_Fit_Model_Peak
end

parameter_size(::Gaussian_Peak) = 3

function eval_y(::Gaussian_Peak,x::Real,θ::AbstractVector)
    h=θ[1]
    μ=θ[2]
    σ=θ[3]

    t = ((x-μ)/σ)^2
    
    h*exp(-t/2)
end

# ----------------------------------------------------------------

using StaticArrays

struct Peak_Motif{P<:Abstract_Fit_Model_Peak} <: Abstract_Fit_Model_Peak
    _peak::P
    _profile::Matrix{Float64} # position, height
end

# -- Specialization for Gaussian Peak
#    θ = (h,σ)
#
parameter_size(pm::Peak_Motif) = 2


# function eval_y(pm::Peak_Motif{Gaussian_Peak},x::Real,θ::AbstractVector{T}) where {T}
#     @assert length(θ) == parameter_size(pm)

#     h_glob = θ[1]
#     σ_glob = θ[2]

#     θ_loc    = Vector{T}(undef,3)
#     θ_loc[3] = σ_glob

#     n::Int = size(pm._profile,1)

#     sum = 0
#     for i in 1:n
#         μ_loc = pm._profile[i,1]
#         h_loc = pm._profile[i,2]

#         θ_loc[1] = h_glob*h_loc
#         θ_loc[2] = μ_loc
        
#         sum += eval_y(pm._peak, x, θ_loc)
#     end

#     sum
# end
function eval_y(pm::Peak_Motif{Gaussian_Peak},x::Real,θ::AbstractVector{T}) where {T}
    @assert length(θ) == parameter_size(pm)

    h_glob = θ[1]
    σ_glob = θ[2]

    n = size(pm._profile,1)

    sum = 0
    for i in 1:n
        μ_loc = pm._profile[i,1]
        h_loc = pm._profile[i,2]

        sum += eval_y(pm._peak, x, @SVector T[ h_glob*h_loc,μ_loc,σ_glob])
    end

    sum
end 
 

# ----------------------------------------------------------------

using ForwardDiff, NLS_Solver, LinearAlgebra, BenchmarkTools

import NLS_Solver


# Using "template" parameters allows to have 1 alloc in eval_r (versus
# 6 if one uses X,Y::AbstractVector)
#
struct NLS_ForwardDiff_From_Fit_Model{X_ELEMENT_TYPE,
                                      Y_ELEMENT_TYPE,
                                      X_TYPE <: AbstractVector{X_ELEMENT_TYPE},
                                      Y_TYPE <: AbstractVector{Y_ELEMENT_TYPE}} <: NLS_Solver.AbstractNLS
    _fit_model::Abstract_Fit_Model
    _X::X_TYPE
    _Y::Y_TYPE
 
    # if X is a n x m 2D array, we assume that each row is an
    # evaluation site (a point of R^m), hence the computed Y is of length n = size(X,1)
    #
    function NLS_ForwardDiff_From_Fit_Model(fit_model::Abstract_Fit_Model,
                                            X::X_TYPE,Y::Y_TYPE) where{X_ELEMENT_TYPE,
                                                                       Y_ELEMENT_TYPE,
                                                                       X_TYPE <: AbstractVector{X_ELEMENT_TYPE},
                                                                       Y_TYPE <: AbstractVector{Y_ELEMENT_TYPE}} 
        @assert length(X) == length(Y)

        new{X_ELEMENT_TYPE,Y_ELEMENT_TYPE,X_TYPE,Y_TYPE}(fit_model,X,Y)
    end
end


NLS_Solver.parameter_size(nls::NLS_ForwardDiff_From_Fit_Model) = parameter_size(nls._fit_model)
NLS_Solver.residue_size(nls::NLS_ForwardDiff_From_Fit_Model)  = length(nls._Y)

function NLS_Solver.eval_r(nls::NLS_ForwardDiff_From_Fit_Model,θ::AbstractVector) 
    map(((X_i,Y_i);)->Y_i-eval_y(nls._fit_model,X_i,θ),zip(nls._X,nls._Y))
end


function NLS_Solver.eval_r_J(nls::NLS_ForwardDiff_From_Fit_Model, θ::AbstractVector{T}) where {T}
    
    r_evaluation = (r,θ)->(r .= eval_r(nls,θ))
    
    r = Vector{T}(undef,residue_size(nls))

    J = ForwardDiff.jacobian(r_evaluation, r, θ)

    r,J
end

# ================================================================
# Demo
# ================================================================

n = 10
X = Float64[1:n;]
Y = rand(n)

model = Gaussian_Peak()

θ = rand(parameter_size(model))

nls = NLS_ForwardDiff_From_Fit_Model(model,X,Y)


eval_r(nls,θ)
eval_r_J(nls,θ)

conf = Levenberg_Marquardt_Conf()
result=solve(nls, θ, conf)

norm(eval_r(nls,θ),2)
norm(eval_r(nls,solution(result)),2)

# ----------------------------------------------------------------

profile=rand(5,2)

model = Peak_Motif(Gaussian_Peak(),profile)

θ = rand(parameter_size(model))

@btime eval_y(model,2.0,θ)
       
nls = NLS_ForwardDiff_From_Fit_Model(model,X,Y)

eval_r(nls,θ)
@btime eval_r($nls,$θ)
@btime eval_r_J(nls,θ)
@btime solve($nls, $θ, $conf)

# result=solve(nls, θ, conf)


