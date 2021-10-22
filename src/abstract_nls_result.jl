export converged, iteration_count, objective_value, solution

# ================================================================
# Result abstraction for unconstrained problem
# ================================================================
#
@doc raw"""
```julia
abstract type AbstractNLSResult end
```

NLS solver result abstraction
"""
abstract type AbstractNLSResult end

"""
    converged(::AbstractNLSResult)

Return `true` if the solver converged
"""
converged(::AbstractNLSResult) = error("To implement")

"""
    iteration_count(::AbstractNLSResult)

Return the number of consumed iteration
"""
iteration_count(::AbstractNLSResult) = error("To implement")


"""
    objective_value(::AbstractNLSResult)

Returns objective value at the point [`solution`](@ref).
"""
objective_value(r::AbstractNLSResult) =  error("To implement")

"""
    solution(::AbstractNLSResult)

Returns the founded solution 
"""
solution(r::AbstractNLSResult) =  error("To implement")


# ================================================================
# Result abstraction for bound constrained problems
# ================================================================
#
# For the moment only use an alias as we do not have specific methods
#
@doc raw"""
```julia
abstract type AbstractNLSResult end
```

NLS BC solver result abstraction
"""
const AbstractNLSBCResult = AbstractNLSResult

