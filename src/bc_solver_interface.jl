export solve

"""
```julia
solve(nls::AbstractNLS,
      θ_init::AbstractVector,
      bc::BoundConstraints,
      conf::Abstract_BC_Solver_Conf) -> Abstract_Solver_Result
```

Generic interface to solve a [`AbstractNLS`](@ref) problem **with
bound constraints**.

The used algorithm is defined through [`Abstract_BC_Solver_Conf`](@ref)
specializations.

The method returns a [`Abstract_BC_Solver_Result`](@ref) specialization.
"""
solve(nls::AbstractNLS,
      θ_init::AbstractVector,
      bc::BoundConstraints,
      conf::Abstract_BC_Solver_Conf) = @assert(false,"To implement")

               
