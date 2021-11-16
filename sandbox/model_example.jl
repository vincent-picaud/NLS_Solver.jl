using ForwardDiff, NLS_Solver, LinearAlgebra, BenchmarkTools, StaticArrays

#
# A model example
#
abstract type Abstract_Model2Fit end

parameter_size(::Abstract_Model2Fit) = @assert(false,"To implement!")

# Here eval Y_i for point X_i. When performing fit we have a collection of such X_i
#
eval_y(::Abstract_Model2Fit,X_i::Any,θ::AbstractVector)  = @assert(false,"To implement!")

# ----------------------------------------------------------------

abstract type Abstract_Model2Fit_Peak <: Abstract_Model2Fit end

# ----------------------------------------------------------------

struct Gaussian_Peak <: Abstract_Model2Fit_Peak
end

parameter_size(::Gaussian_Peak) = 3

function eval_y(m::Gaussian_Peak,x::Real,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)
    
    h=θ[1]
    μ=θ[2]
    σ=θ[3]

    t = ((x-μ)/σ)^2
    
    h*exp(-t/2)
end

# ----------------------------------------------------------------

#
# julia> @btime eval_y($model,2.0,$θ)
#   281.779 ns (8 allocations: 208 bytes)
# 0.06722452348922268
# 
# struct Model2Fit_Sum <: Abstract_Model2Fit
#     _left::Abstract_Model2Fit
#     _right::Abstract_Model2Fit
# end
#
# versus:
#
# struct Model2Fit_Sum{LEFT_TYPE<:Abstract_Model2Fit,
#                      RIGHT_TYPE<:Abstract_Model2Fit} <: Abstract_Model2Fit
#     _left::LEFT_TYPE
#     _right::RIGHT_TYPE
# end
#
# julia> @btime eval_y($model,2.0,$θ)
#   77.596 ns (0 allocations: 0 bytes)
# 0.018749729970880733
#
struct Model2Fit_Sum{LEFT_TYPE<:Abstract_Model2Fit,
                     RIGHT_TYPE<:Abstract_Model2Fit} <: Abstract_Model2Fit
    _left::LEFT_TYPE
    _right::RIGHT_TYPE
end

import Base: (+)

Base. +(left::Abstract_Model2Fit,right::Abstract_Model2Fit) = Model2Fit_Sum(left,right)

parameter_size(m::Model2Fit_Sum) = parameter_size(m._left) + parameter_size(m._right)

function eval_y(m::Model2Fit_Sum,x::Any,θ::AbstractVector)
    @assert length(θ) == parameter_size(m)
    
    eval_y(m._left,x,@view θ[1:parameter_size(m._left)]) +
        eval_y(m._right,x,@view θ[(parameter_size(m._left)+1):end])
end

# test
model = Gaussian_Peak() + Gaussian_Peak()
θ=rand(parameter_size(model))
@btime eval_y($model,2.0,$θ)

# ----------------------------------------------------------------



struct Peak_Motif{P<:Abstract_Model2Fit_Peak} <: Abstract_Model2Fit_Peak
    _peak::P
    _profile::Matrix{Float64} # position, height
end

# -- Specialization for Gaussian Peak
#    θ = (h,σ)
#
parameter_size(pm::Peak_Motif) = 2


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



import NLS_Solver


# Important: using "template" parameters allows to have 1 alloc in
# eval_r (versus 6 if one uses X,Y::AbstractVector)
#
#
# Finally MODEL2FIT_TYPE <: Abstract_Model2Fit is parametrized 
# julia> @btime eval_r($nls,$θ)
#   2.724 μs (44 allocations: 928 bytes)
# 10-element Vector{Float64}:
#  0.8056529547193196
#  0.10344414895180448
#  0.7746398216971467
#  0.5705008291976804
#  0.758504377425209
#  0.6572756090278895
#  0.22021418677443244
#  0.9012363674393303
#  0.25785439137881183
#  0.7594027591650859
#
# julia> @btime eval_r($nls,$θ)
#   842.590 ns (1 allocation: 160 bytes)
# 10-element Vector{Float64}:
#  -0.6484862692515055
#   0.10436137221132202
#   0.3416337767340343
#  -0.5195066030780292
#   0.0635334492635572
#  -0.2671004288239941
#  -0.13064413518496942
#  -0.143200160754688
#   0.517014266738562
#   0.21068675449678864
# 
struct NLS_ForwardDiff_From_Fit_Model{MODEL2FIT_TYPE <: Abstract_Model2Fit,
                                      X_ELEMENT_TYPE,
                                      Y_ELEMENT_TYPE,
                                      X_TYPE <: AbstractVector{X_ELEMENT_TYPE},
                                      Y_TYPE <: AbstractVector{Y_ELEMENT_TYPE}} <: NLS_Solver.AbstractNLS
    _fit_model::MODEL2FIT_TYPE
    _X::X_TYPE
    _Y::Y_TYPE
 
    function NLS_ForwardDiff_From_Fit_Model(fit_model::MODEL2FIT_TYPE,
                                            X::X_TYPE,Y::Y_TYPE) where{MODEL2FIT_TYPE <: Abstract_Model2Fit,
                                                                       X_ELEMENT_TYPE,
                                                                       Y_ELEMENT_TYPE,
                                                                       X_TYPE <: AbstractVector{X_ELEMENT_TYPE},
                                                                       Y_TYPE <: AbstractVector{Y_ELEMENT_TYPE}} 
        @assert length(X) == length(Y)

        new{MODEL2FIT_TYPE,X_ELEMENT_TYPE,Y_ELEMENT_TYPE,X_TYPE,Y_TYPE}(fit_model,X,Y)
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# An attempt of wrapping using Zygote
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Does not work at all in our case:
# julia> @btime eval_r_J($nls_zygote,$θ)
#   1.421 s (2632107 allocations: 81.61 MiB) 
#
# versus Forward diff
#
# julia> @btime eval_r_J($nls,$θ)
#   36.179 μs (9 allocations: 26.11 KiB)
#
# using Zygote

# struct NLS_Zygote_From_Fit_Model{MODEL2FIT_TYPE <: Abstract_Model2Fit,
#                                  X_ELEMENT_TYPE,
#                                  Y_ELEMENT_TYPE,
#                                  X_TYPE <: AbstractVector{X_ELEMENT_TYPE},
#                                  Y_TYPE <: AbstractVector{Y_ELEMENT_TYPE}} <: NLS_Solver.AbstractNLS
#     _fit_model::MODEL2FIT_TYPE
#     _X::X_TYPE
#     _Y::Y_TYPE
 
#     function NLS_Zygote_From_Fit_Model(fit_model::MODEL2FIT_TYPE,
#                                             X::X_TYPE,Y::Y_TYPE) where{MODEL2FIT_TYPE <: Abstract_Model2Fit,
#                                                                        X_ELEMENT_TYPE,
#                                                                        Y_ELEMENT_TYPE,
#                                                                        X_TYPE <: AbstractVector{X_ELEMENT_TYPE},
#                                                                        Y_TYPE <: AbstractVector{Y_ELEMENT_TYPE}} 
#         @assert length(X) == length(Y)

#         new{MODEL2FIT_TYPE,X_ELEMENT_TYPE,Y_ELEMENT_TYPE,X_TYPE,Y_TYPE}(fit_model,X,Y)
#     end
# end


# NLS_Solver.parameter_size(nls::NLS_Zygote_From_Fit_Model) = parameter_size(nls._fit_model)
# NLS_Solver.residue_size(nls::NLS_Zygote_From_Fit_Model)  = length(nls._Y)

# function NLS_Solver.eval_r(nls::NLS_Zygote_From_Fit_Model,θ::AbstractVector{T}) where T
#     # Cannot use adjoint as zip() currently not supported 
#     #    map(((X_i,Y_i);)->Y_i-eval_y(nls._fit_model,X_i,θ),zip(nls._X,nls._Y))
#     r = map(X_i->eval_y(nls._fit_model,X_i,θ),nls._X)
#     nls._Y-r
# end


# function NLS_Solver.eval_r_J(nls::NLS_Zygote_From_Fit_Model, θ::AbstractVector{T}) where {T}
    
#     r = eval_r(nls,θ)

#     J = jacobian(θ->eval_r(nls,θ), θ)

#     r,J
# end

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# An attempt of wrapping using Enzyme
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Cannot make it work...
#
# julia> autodiff(t->eval_y(model,2.0,t),θ)
# ┌ Warning: ("reverse differentiating jl_apply_generic call without split mode", eval_y, 0x00000003)
# └ @ Enzyme.Compiler ~/.julia/packages/Enzyme/afnXq/src/compiler.jl:260
# ┌ Warning: ("done reverse differentiating jl_apply_generic call without split mode", eval_y, ())
# └ @ Enzyme.Compiler ~/.julia/packages/Enzyme/afnXq/src/compiler.jl:325
# ()


# ================================================================
# Demo
# ================================================================

X = rand(100)
Y = rand(100)

model = Gaussian_Peak()+Gaussian_Peak()+Gaussian_Peak()

θ = rand(parameter_size(model))

nls = NLS_ForwardDiff_From_Fit_Model(model,X,Y)

@btime eval_r($nls,$θ)
@btime eval_r_J($nls,$θ)

# ----------------------------------------------------------------

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


