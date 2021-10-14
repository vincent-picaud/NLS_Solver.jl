export BoundConstraints
export lower_bound,upper_bound
export project!

using ReadOnlyArrays

"""
Store bound constraints ``[l,u]``

Presence of `NaN` component and the ``l\\le u`` condition is checked
at construction time. Note however that some components can be
infinite.

```jldoctest
julia> bc = BoundConstraints(5)
BoundConstraints{Float64, 1, Vector{Float64}, Vector{Float64}}([0.0, 0.0, 0.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0, 1.0])

```
"""
struct BoundConstraints{ELT<:Real,N,LBT<:AbstractArray{ELT,N},UBT<:AbstractArray{ELT,N}}
    _lb::LBT
    _ub::UBT

    function BoundConstraints(lb::AbstractArray{ELT,N},ub::AbstractArray{ELT,N}) where{ELT,N}
        @assert size(lb)==size(ub)
        # note: also fails in case of NaN
        @assert all((lb_i ≤ ub_i for (lb_i,ub_i) ∈ zip(lb,ub)))
        new{ELT,N,typeof(lb),typeof(ub)}(lb,ub)
    end
end 

BoundConstraints(n::Int) = BoundConstraints(zeros(n),ones(n))
BoundConstraints(T::Type,n::Int) = BoundConstraints(zeros(T,n),ones(T,n))


import Base: size, length, axes, in
"""
    axes(bc::BoundConstraints)

Return bound axes
"""
Base.axes(bc::BoundConstraints) = axes(bc._lb)

"""
    length(bc::BoundConstraints)

Return bound length
"""
Base.length(bc::BoundConstraints) = length(bc._lb)

"""
    size(bc::BoundConstraints)

Return bound size
"""
Base.size(bc::BoundConstraints) = size(bc._lb)

"""
    in(bc::BoundConstraints)

Check if ``x\\in [l,u]``
"""
function in(x::AbstractArray{<:Real,N},bc::BoundConstraints{<:Real,N}) where {N}
    size(x)==size(bc) && all(zip(bc._lb,x,bc._ub)) do (lbi,xi,ubi) begin lbi ≤ xi ≤ ubi end end
end 

"""
    lower_bound(bc::BoundConstraints)

Return lower bound `l`
"""
lower_bound(bc::BoundConstraints) = ReadOnlyArray(bc._lb)

"""
    upper_bound(bc::BoundConstraints)

Return upper bound `u`
"""
upper_bound(bc::BoundConstraints) = ReadOnlyArray(bc._ub)

"""
    project!(x::AbstractArray{<:Real,N},bc::BoundConstraints{<:Real,N})

Project `x` such that ``x \\in [l,u]`` is fullfiled.

"""
function project!(x::AbstractArray{<:Real,N},bc::BoundConstraints{<:Real,N}) where {N}
    @assert size(x) == size(bc)
    for (idx,(lb_idx,ub_idx)) in enumerate(zip(bc._lb,bc._ub))
        @inbounds x[idx]=max(lb_idx,min(x[idx],ub_idx))
    end
    x
end 