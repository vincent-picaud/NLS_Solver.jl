export AbstractQuadSolverResult
export converged, iteration_count,objective_value,multiplier_τ,solution

"""
    abstract type AbstractQuadSolverResult end

Quadratic solver result abstraction

"""
abstract type AbstractQuadSolverResult end


"""
    converged(::AbstractQuadSolverResult)

Return `true` if the solver converged
"""
converged(::AbstractQuadSolverResult) = error("To implement")

"""
    iteration_count(::AbstractQuadSolverResult)

Return the number of consumed iteration
"""
iteration_count(::AbstractQuadSolverResult) = error("To implement")


"""
    objective_value(::AbstractQuadSolverResult)

Returns objective value at the point [`solution`](@ref).
"""
objective_value(r::AbstractQuadSolverResult) = r._fobj # TODO: remove this default! Tue Oct 19 10:30:46 AM CEST 2021

"""
    multiplier_τ(::AbstractQuadSolverResult)

Returns the multipliers stored in a compact form (see τ definition, TODO)
"""
multiplier_τ(r::AbstractQuadSolverResult) = ReadOnlyArray(r._τ) # TODO: remove this default! Tue Oct 19 10:30:46 AM CEST 2021

"""
    solution(::AbstractQuadSolverResult)

Returns the founded solution 
"""
solution(r::AbstractQuadSolverResult) = ReadOnlyArray(r._x) # TODO: remove this default! Tue Oct 19 10:30:46 AM CEST 2021
