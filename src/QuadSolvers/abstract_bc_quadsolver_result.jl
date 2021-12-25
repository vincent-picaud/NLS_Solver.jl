"""
    abstract type Abstract_BC_QuadSolver_Result end

Quadratic solver result abstraction

"""
abstract type Abstract_BC_QuadSolver_Result end


"""
    converged(::Abstract_BC_QuadSolver_Result)

Return `true` if the solver converged
"""
converged(::Abstract_BC_QuadSolver_Result) = @assert(false,"To implement")

"""
    iteration_count(::Abstract_BC_QuadSolver_Result)

Return the number of consumed iteration
"""
iteration_count(::Abstract_BC_QuadSolver_Result) = @assert(false,"To implement")


"""
    objective_value(::Abstract_BC_QuadSolver_Result)

Returns objective value at the point [`solution`](@ref).
"""
objective_value(r::Abstract_BC_QuadSolver_Result) = @assert(false,"To implement")

"""
    multiplier_τ(::Abstract_BC_QuadSolver_Result)

Returns the multipliers stored in a compact form (see τ definition, TODO)
"""
multiplier_τ(r::Abstract_BC_QuadSolver_Result) = @assert(false,"To implement")

"""
    solution(::Abstract_BC_QuadSolver_Result)

Returns the founded solution 
"""
solution(r::Abstract_BC_QuadSolver_Result) = @assert(false,"To implement")
