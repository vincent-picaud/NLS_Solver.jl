export solve

"""
```julia
solve(nls::AbstractNLS,
      θ_init::AbstractVector{<:Real},
      bc::BoundConstraints{<:Real},
      conf::AbstractNLSBCConf) -> AbstractNLSResult
```

Generic interface to solve a [`AbstractNLS`](@ref) problem **with
bound constraints**.

The algorithm to be used is defined through `conf` of type
[`AbstractNLSBCConf`](@ref) specializations.

The method returns a [`AbstractNLSResult`](@ref) specialization.
"""
solve(nls::AbstractNLS,
      θ_init::AbstractVector,
      bc::BoundConstraints{<:Real},
      conf::AbstractNLSBCConf) = error("To implement")

               
