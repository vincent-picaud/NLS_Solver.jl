# Some misc functions related to quadratic optimization problems
#
import LinearAlgebra: Symmetric

@doc raw"""
```julia
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

Using the previous results, this function returns:
```math
\max \mid x^\star - P_{[l,u]}(x^\star-(Q.x^\star+q)) \mid
```
for a local optimal quantity this must be null 
"""
function check_first_order(Q::Symmetric{<:Real},
                           q::AbstractVector{<:Real},
                           xstar::AbstractVector{<:Real},
                           bc::BoundConstraints{<:Real,1})
    v = xstar-(Q*xstar+q)
    v = project!(v,bc)

    maximum(abs.(xstar .- v))
end 
