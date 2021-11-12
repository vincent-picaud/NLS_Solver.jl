export solve

"""
```julia
solve(nls::AbstractNLS,
      θ_init::AbstractVector{<:Real},
      conf::AbstractNLSConf) -> AbstractNLSResult
```

Generic interface to solve a [`AbstractNLS`](@ref) problem.

The algorithm to be used is defined through `conf` of type
[`AbstractNLSConf`](@ref) specializations.

The method returns a [`AbstractNLSResult`](@ref) specialization.
"""
solve(nls::AbstractNLS,
      θ_init::AbstractVector,
      conf::AbstractNLSConf) = @assert(false,"To implement")

               
