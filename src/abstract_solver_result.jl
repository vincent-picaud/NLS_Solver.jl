export Abstract_Solver_Result
export converged, iteration_count, objective_value, solution

# Result abstraction for unconstrained problem ================
#
@doc raw"""
```julia
abstract type Abstract_Solver_Result end
```

This is the base type returned by the [`solve`](@ref) method. It
contains the information related to the founded solution.

# Interface
- [`converged`](@ref) 
- [`iteration_count`](@ref) 
- [`objective_value`](@ref) 
- [`solution`](@ref) 

# Implementations
- [`LevenbergMarquardt_Result`](@ref) 
- [`LevenbergMarquardt_BC_Result`](@ref) 

"""
abstract type Abstract_Solver_Result end

"""
    converged(::Abstract_Solver_Result)

Return `true` if the solver converged

See: [`Abstract_Solver_Result`](@ref) 

"""
converged(::Abstract_Solver_Result) = @assert(false,"To implement")

"""
    iteration_count(::Abstract_Solver_Result)

Return the number of consumed iteration

See: [`Abstract_Solver_Result`](@ref) 

"""
iteration_count(::Abstract_Solver_Result) = @assert(false,"To implement")


"""
    objective_value(::Abstract_Solver_Result)

Returns objective value at the point [`solution`](@ref).

See: [`Abstract_Solver_Result`](@ref) 

"""
objective_value(r::Abstract_Solver_Result) =  @assert(false,"To implement")

"""
    solution(::Abstract_Solver_Result)

Returns the founded solution 

See: [`Abstract_Solver_Result`](@ref) 

"""
solution(r::Abstract_Solver_Result) =  @assert(false,"To implement")


