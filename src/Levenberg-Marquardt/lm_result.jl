@doc raw"""

The structure returned by [`solve`](@ref) when using the
[`LevenbergMarquardt_Conf`](@ref) method.

See [`Abstract_Solver_Result`](@ref) 
"""
struct LevenbergMarquardt_Result{T<:Real} <: Abstract_Solver_Result
    _converged::Bool
    _iter_count::Int
    _fobj::T
    _solution::Vector{T}

    # Do not use Base.@kwdef for Julia 1 compatibility reason
    function LevenbergMarquardt_Result(;
                                       _converged,
                                       _iter_count,
                                       _fobj,
                                       _solution::Vector{T}) where {T<:Real}

        new{T}(_converged, _iter_count, _fobj, _solution)
    end

end 
converged(r::LevenbergMarquardt_Result) = r._converged
iteration_count(r::LevenbergMarquardt_Result) = r._iter_count
objective_value(r::LevenbergMarquardt_Result) = r._fobj
solution(r::LevenbergMarquardt_Result) = ReadOnlyArray(r._solution)
