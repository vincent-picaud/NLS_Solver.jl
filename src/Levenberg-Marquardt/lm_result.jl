# The structure to be returned, used in: lm.jl
# Also see: lm_conf.jl
#
Base.@kwdef struct LevenbergMarquardt_Result{T<:Real} <: AbstractNLSResult
    _converged::Bool
    _iter_count::Int
    _fobj::T
    _solution::Vector{T}
end 
converged(r::LevenbergMarquardt_Result) = r._converged
iteration_count(r::LevenbergMarquardt_Result) = r._iter_count
objective_value(r::LevenbergMarquardt_Result) = r._fobj
solution(r::LevenbergMarquardt_Result) = ReadOnlyArray(r._solution)
