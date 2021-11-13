
@doc raw"""
Compute δL = L(0)-L(h) where L is the quadratic model

```math
L(h)=\frac{1}{2}\| r(θ) \|_2^2 + \langle \nabla f, h \rangle + \frac{1}{2}\langle \nabla^2 f h, h \rangle
```

with ``f(\theta)=\frac{1}{2}\| r(θ) \|_2^2``, ``\nabla f = J^t r`` and ``\nabla^2 f = J^t J``

A direct computation gives:


```math
δL = L(0)-L(h) = -\left(  \langle J^tr, h \rangle + \frac{1}{2}\langle \nabla^2 f h, h \rangle \right)
```

However one can avoid the computation of ``\nabla^2 f h`` if one uses the
fact that ``h`` is solution of:

```math
(\nabla^2 f + \mu I)h + \nabla f = 0
```

With this hypothesis, one gets:

```math
δL = L(0)-L(h) = \frac{1}{2} \langle h, \mu h - \nabla f \rangle
```
"""
function compute_δL_unconstrained(∇f::AbstractVector,
                                  μ::Real,
                                  h::AbstractVector)

    # The syntax to desctructure is unexpected... 
    # https://discourse.julialang.org/t/argument-destructuring-and-anonymous-functions/24893/2
    #
    mapreduce(((h_i,∇f_i),)->h_i*(h_i-μ*∇f_i),+,zip(h,∇f))/2
end


@doc raw"""

Same idea than [`compute_δL_unconstrained`](@ref), however when bound
constraints are present ``h`` is such that:

```math
(\nabla^2 f + \mu I)h + \nabla f + \tau = 0
```
it follows that:

```math
δL = L(0)-L(h) = \frac{1}{2} \langle h, \mu h + \tau - \nabla f \rangle
```
"""
function compute_δL_constrained(∇f::AbstractVector,
                                μ::Real,
                                τ::AbstractVector,
                                h::AbstractVector)
    mapreduce(((h_i,τ_i,∇f_i),)->h_i*(h_i+τ_i-μ*∇f_i),+,zip(h,τ,∇f))/2
end


@doc raw"""

Compute true variation of the real model: ``δf = \frac{1}{2}(r^t(θ)r(θ)-r^t(θ+h)r(θ+h))``

Contrary to ``δL`` things are simpler. However a trick is to use an equivalent formulation:
```math
δf = \frac{1}{2}(r^t(θ)r(θ)-r^t(θ+h)r(θ+h)) = \frac{1}{2}(r(θ)-r(θ+h))^t(r(θ)+r(θ+h))
```
that has a better numerical behavior.
"""
function compute_δf(r::AbstractVector,
                    r_new::AbstractVector)
    mapreduce(((r_i,r_new_i),)->(r_i-r_new_i)*(r_i+r_new_i),+,zip(r,r_new))
end
