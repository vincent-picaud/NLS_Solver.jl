export converged, iteration_count, objective_value, solution

# ================================================================
# Result abstraction for unconstrained problem
# ================================================================
#
@doc raw"""
```julia
abstract type Abstract_Solver_Result end
```

NLS solver result abstraction
"""
abstract type Abstract_Solver_Result end

"""
    converged(::Abstract_Solver_Result)

Return `true` if the solver converged
"""
converged(::Abstract_Solver_Result) = @assert(false,"To implement")

"""
    iteration_count(::Abstract_Solver_Result)

Return the number of consumed iteration
"""
iteration_count(::Abstract_Solver_Result) = @assert(false,"To implement")


"""
    objective_value(::Abstract_Solver_Result)

Returns objective value at the point [`solution`](@ref).
"""
objective_value(r::Abstract_Solver_Result) =  @assert(false,"To implement")

"""
    solution(::Abstract_Solver_Result)

Returns the founded solution 
"""
solution(r::Abstract_Solver_Result) =  @assert(false,"To implement")


