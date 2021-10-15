export AbstractRegularizationSchedule
export regularization_factor, burning_phase

export NoRegularizationSchedule, ExpRegularizationSchedule

""" 
    abstract type AbstractRegularizationSchedule end

Generate a finite sequence of reals that are getting closer and closer
to Float64(1).

These reals are interpreted as regularization factors, at the end
(when [`burning_phase`](@ref)=false there is no more regularization
([`regularization_factor`](@ref) = 1).
"""
abstract type AbstractRegularizationSchedule end

"""
    regularization_factor(rs::AbstractRegularizationSchedule, iter::Int) 

Return the regularization factor value
"""
regularization_factor(rs::AbstractRegularizationSchedule, iter::Int) = error("to implement")


"""
    burning_phase(rs::AbstractRegularizationSchedule, iter::Int)

Return a Boolean indicating if we are still in the burning phase.
""" 
burning_phase(rs::AbstractRegularizationSchedule, iter::Int) = regularization_factor(rs,iter)!=Float64(1)

# ================================================================
# Some concrete types...
# ================================================================

"""
```julia
struct NoRegularizationSchedule <: AbstractRegularizationSchedule
  ...
end
```
No regularization. The returned factor is always the unity and the
`burning_phase` predicate is always false.
"""
struct NoRegularizationSchedule <: AbstractRegularizationSchedule
end

regularization_factor(rs::NoRegularizationSchedule, iter::Int) = Float64(1)

# ================================================================

@doc raw"""
```julia
struct ExpRegularizationSchedule <: AbstractRegularizationSchedule
  ...
end
```

Exponentially decreasing factor at iteration `k` is 
```math
\alpha(k) = \left\{\begin{\array}{ll}
c \exp{-\log{c}\frac{k-1}{n-1}}, & \text{if } k<n \\
1, & \text{if } k \ge n \\
\end{array}\right.
```
By construction we have `α(1)=c` and `α(k≥n)=1`.
"""
struct ExpRegularizationSchedule <: AbstractRegularizationSchedule
    c::Float64
    n::Int
    _β::Float64

    function ExpRegularizationSchedule(c::Real,n::Integer)
        @assert c≥1
        @assert n>1

        # β is such that: α(k) = c \exp{β*(k-1)}
        _β = -log(c)/(n-1)

        new(c,n,_β)
    end
end

burning_phase(rs::ExpRegularizationSchedule, iter::Int) = iter ≤ rs.n
regularization_factor(rs::ExpRegularizationSchedule, iter::Int) = rs.c*exp(rs._β*(iter-1))
