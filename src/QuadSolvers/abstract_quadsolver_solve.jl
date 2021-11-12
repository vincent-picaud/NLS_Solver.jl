export solve

using LinearAlgebra: Symmetric


@doc raw"""
```julia
solve(Q::Symmetric{<:Real},
      q::AbstractVector{<:Real},
      x_init::AbstractVector{<:Real},
      bc::BoundConstraints{<:Real,1},
      conf::AbstractQuadSolverConf)::AbstractQuadSolverResult
```


Generic interface to solve a quadratic optimization problem of the form:

```math
\min\limits_{\mathbf{a}\le \mathbf{x} \le \mathbf{b}} \frac{1}{2} \mathbf{x}^t\mathbf{Q}\mathbf{x} + \mathbf{q}^t\mathbf{x}
```

The algorithm to be used is defined through
[`AbstractQuadSolverConf`](@ref) specializations. The method returns a
[`AbstractQuadSolverResult`](@ref) specialization.

"""
solve(Q::Symmetric{<:Real},
      q::AbstractVector{<:Real},
      x_init::AbstractVector{<:Real},
      bc::BoundConstraints{<:Real,1},
      conf::AbstractQuadSolverConf) = @assert(false,"To implement")

               
