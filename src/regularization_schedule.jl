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

Exponentially decreasing factor. 

At iteration `k` is the factor α is:
```math
\alpha(k) = \left\{\begin{array}{ll}
c \exp{\left(-\log{c}\left(\frac{k-1}{n}\right)\right)}, & \text{if } k \le n \\
1, & \text{if } k \ge n \\
\end{array}\right.
```
where n is the last burning phase iteration (id the last iteration for which α>1).

By construction we have α(1)=c and α(k>n)=1.
"""
struct ExpRegularizationSchedule <: AbstractRegularizationSchedule
    factor::Float64
    burning_last_iter::Int
    _β::Float64

    function ExpRegularizationSchedule(;factor::Real,burning_last_iter::Integer)
        @assert factor ≥ 1
        @assert burning_last_iter>1

        # β is such that: α(k) = c \exp{β*k}
        _β = -log(factor)/burning_last_iter

        new(factor,burning_last_iter,_β)
    end
end

function regularization_factor(rs::ExpRegularizationSchedule, iter::Int)
    ifelse(iter <=  rs.burning_last_iter, rs.factor*exp(rs._β*(iter-1)), Float64(1))
end


