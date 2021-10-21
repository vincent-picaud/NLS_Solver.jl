@doc raw"""

A structure to store:
- [`Levenberg_Marquardt`](@ref)
- [`Levenberg_Marquardt_BC`](@ref)
results.

"""
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
