export solve

"""
```julia
solve(nls::AbstractNLS,
      θ_init::AbstractVector{<:Real},
      conf::Abstract_Solver_Conf) -> Abstract_Solver_Result
```

Generic interface to solve a [`AbstractNLS`](@ref) problem.

The algorithm to be used is defined through `conf` of type
[`Abstract_Solver_Conf`](@ref) specializations.

The method returns a [`Abstract_Solver_Result`](@ref) specialization.
"""
solve(nls::AbstractNLS,
      θ_init::AbstractVector,
      conf::Abstract_Solver_Conf) = @assert(false,"To implement")

               
