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
converged(::AbstractNLSResult) = @assert(false,"To implement")

"""
    iteration_count(::AbstractNLSResult)

Return the number of consumed iteration
"""
iteration_count(::AbstractNLSResult) = @assert(false,"To implement")


"""
    objective_value(::AbstractNLSResult)

Returns objective value at the point [`solution`](@ref).
"""
objective_value(r::AbstractNLSResult) =  @assert(false,"To implement")

"""
    solution(::AbstractNLSResult)

Returns the founded solution 
"""
solution(r::AbstractNLSResult) =  @assert(false,"To implement")


