export AbstractQuadSolverResult
export converged, objective_value, multiplier_τ, iteration_count, solution
    
"""

Abstract away the result of the quadratic solver

TODO example

# Interface

"""
abstract type AbstractQuadSolverResult end


"""
    converged(::AbstractQuadSolverResult)

Returns `true` if the solver converged
"""
converged(::AbstractQuadSolverResult) = error("To implement")

"""
    iteration_count(::AbstractQuadSolverResult)

Returns the number of iterations
"""
iteration_count(::AbstractQuadSolverResult) = error("To implement")

"""
    objective_value(::AbstractQuadSolverResult)

Returns the final objectif function value
"""
objective_value(::AbstractQuadSolverResult) = error("To implement")

"""
    multiplier_τ(::AbstractQuadSolverResult)

Returns the Lagrange multipliers encoded in the `τ` array:
- `τ[i]=0` the constraint is inactive
- `τ[i]<0` the lower bound constraint is active, its mulitplier is `λ[i] = -τ[i]`
- `τ[i]>0` the upper bound constraint is active, its mulitplier is `μ[i] = +τ[i]`

Also note that, if the algorithm converged, then `∇objf = -τ` 
"""
multiplier_τ(::AbstractQuadSolverResult) = error("To implement")


"""
    solution(::AbstractQuadSolverResult)

Returns the founded solution 
"""
solution(::AbstractQuadSolverResult) = error("To implement")

