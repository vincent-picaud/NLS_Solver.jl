# Some misc functions related to quadratic optimization problems
#
import LinearAlgebra: Symmetric

@doc raw"""
```julia
check_first_order(∇f::AbstractVector{<:Real},
                  xstar::AbstractVector{<:Real},
                  bc::BoundConstraints{<:Real,1})

check_first_order(Q::Symmetric{<:Real},
                  q::AbstractVector{<:Real},
                  xstar::AbstractVector{<:Real},
                  bc::BoundConstraints{<:Real,1})
```

Check First-Order Conditions 
(see [Bound Constrained Optimization slides](https://wiki.mcs.anl.gov/leyffer/images/0/01/07-bndCons.pdf))

If ``x^\star=\arg\min f(x), x\in[l,u]`` then:

```math
\partial_i f(x^\star) = \left\{\begin{array}{ll}
\ge 0, & \text{if } x^\star[i] = l[i] \\
= 0, & \text{if } l[i] \le x^\star[i] \le u[i] \\
\le 0, & \text{if } x^\star[i] = u[i] \\
\end{array}
\right.
```

This is equivalent to:
```math
x^\star = P_{[l,u]}(x^\star-\nabla f(x^\star))
```

According to the previous result, this function returns:
```math
\max \mid x^\star - P_{[l,u]}(x^\star-(Q.x^\star+q)) \mid
```

For a local stationary point this quantity must be null 

The second function is a wrapper that computes ``∇f=Q.x^\star+q``
"""
function check_first_order(∇f::AbstractVector{<:Real},
                           xstar::AbstractVector{<:Real},
                           bc::BoundConstraints{<:Real,1})

    # the condition assumes feasible constraints
    @assert xstar ∈ bc

    @assert size(∇f)==size(xstar)

    # TODO: refactoring to avoid mem alloc
    # (-> inline computation of max(...sequence...)
    v = xstar-∇f
    v = project!(v,bc)

    maximum(abs.(xstar .- v))
end

function check_first_order(Q::Symmetric{<:Real},
                           q::AbstractVector{<:Real},
                           xstar::AbstractVector{<:Real},
                           bc::BoundConstraints{<:Real,1})
    # TODO: for BlasFloat use BLAS... to avoid unnecessary memory
    # allocs...
    ∇f=Q*xstar+q
    check_first_order(∇f,xstar,bc)
end 
