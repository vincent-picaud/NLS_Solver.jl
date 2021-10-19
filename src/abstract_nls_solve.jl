export solve

"""
```julia
solve(conf::AbstractNLSConf,
      θ_init::AbstractVector{<:Real}) -> AbstractNLSResult
```

Generic interface to solve a [`AbstractNLS`](@ref] problem

The algorithm to be used is defined through
[`AbstractNLSConf`](@ref) specializations. The method returns a
[`AbstractNLSResult`](@ref) specialization.

"""
solve(conf::AbstractNLSConf, θ_init::AbstractVector{<:Real}) = error("To implement")

               
