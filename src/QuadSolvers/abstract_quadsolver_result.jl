export AbstractQuadSolverResult

"""

Abstract away the result of the quadratic solver

TODO example

# Interface

"""
abstract type AbstractQuadSolverResult end


"""

Returns `true` if the solver converged
"""
converged(::AbstractQuadSolverResult) = error("To implement")

"""

Returns the founded solution 
"""
parameters(::AbstractQuadSolverResult) = error("To implement")

