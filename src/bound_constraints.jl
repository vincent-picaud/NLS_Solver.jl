export BoundConstraints
export lower_bound,upper_bound
export project!

using ReadOnlyArrays

"""
Store bound constraints ``[l,u]``

Presence of `NaN` component and the ``l\\le u`` condition is checked
at construction time. Note however that some components can be
infinite.

The following constructors are available:

- Construct ``[0.0,1.0]^n`` 
```julia
BoundConstraints(n)
```
- Construct ``[T(0),T(1)]^n`` where components are of type `T`
```julia
BoundConstraints(T,n)
```
- Construct ``[l,u]`` where `l` and `u` are lower and upper bound vectors
```julia
BoundConstraints(T,n)
```

"""
struct BoundConstraints{ELT<:Real,N,LBT<:AbstractArray{ELT,N},UBT<:AbstractArray{ELT,N}}
    _lb::LBT
    _ub::UBT

    function BoundConstraints(lb::AbstractArray{ELT,N},ub::AbstractArray{ELT,N}) where{ELT,N}
        @assert size(lb)==size(ub)
        @assert axes(lb)==axes(ub)
        # note: also fails in case of NaN
        @assert all((lb_i ≤ ub_i for (lb_i,ub_i) ∈ zip(lb,ub)))
        new{ELT,N,typeof(lb),typeof(ub)}(lb,ub)
    end
end 

BoundConstraints(n::Int) = BoundConstraints(zeros(n),ones(n))
BoundConstraints(T::Type,n::Int) = BoundConstraints(zeros(T,n),ones(T,n))


import Base: eltype, size, length, axes, in


"""
    eltype(bc::BoundConstraints)

Return bound element type

See: [` BoundConstraints`](@ref) 
"""
Base.eltype(bc::BoundConstraints{ELT}) where ELT = ELT


"""
    axes(bc::BoundConstraints)

Return bound axes

See: [` BoundConstraints`](@ref) 
"""
Base.axes(bc::BoundConstraints) = axes(bc._lb) 

"""
    length(bc::BoundConstraints)

Return bound length

See: [` BoundConstraints`](@ref) 
"""
Base.length(bc::BoundConstraints) = length(bc._lb)

"""
    size(bc::BoundConstraints)

Return bound size

See: [` BoundConstraints`](@ref) 
"""
Base.size(bc::BoundConstraints) = size(bc._lb)

"""
    in(bc::BoundConstraints)

Check if ``x\\in [l,u]``

See: [` BoundConstraints`](@ref) 
"""
function in(x::AbstractArray{<:Real,N},bc::BoundConstraints{<:Real,N}) where {N}
    size(x)==size(bc) && all(zip(bc._lb,x,bc._ub)) do (lbi,xi,ubi) begin lbi ≤ xi ≤ ubi end end
end 

"""
    lower_bound(bc::BoundConstraints)

Return lower bound `l`

See: [` BoundConstraints`](@ref) 
"""
lower_bound(bc::BoundConstraints) = ReadOnlyArray(bc._lb)

"""
    upper_bound(bc::BoundConstraints)

Return upper bound `u`

See: [` BoundConstraints`](@ref) 
"""
upper_bound(bc::BoundConstraints) = ReadOnlyArray(bc._ub)

"""
    project!(x::AbstractArray{<:Real,N},bc::BoundConstraints{<:Real,N})

Project `x` such that ``x \\in [l,u]`` is fullfiled.

See: [` BoundConstraints`](@ref) 
"""
function project!(x::AbstractArray{<:Real,N},bc::BoundConstraints{<:Real,N}) where {N}
    @assert size(x) == size(bc)
    for (idx,(lb_idx,ub_idx)) in enumerate(zip(bc._lb,bc._ub))
        @inbounds x[idx]=max(lb_idx,min(x[idx],ub_idx))
    end
    x
end 

# ================================================================

# TODO: define broadcasting interface allowing op like (bc .+ scalar),
# for the moement we have (bc + vector).
# https://docs.julialang.org/en/v1/manual/interfaces/#man-interfaces-broadcasting

import Base: +, -
@doc raw""" 

    Translate bound constraints

```math
[a-τ,b-τ] = [a,b]-τ
```

See: [` BoundConstraints`](@ref) 
"""
Base.:-(bc::BoundConstraints,τ::AbstractArray) = BoundConstraints(lower_bound(bc).-τ,upper_bound(bc).-τ)
    
@doc raw""" 

    Translate bound constraints

```math
[a+τ,b+τ] = [a,b]+τ
```

See: [` BoundConstraints`](@ref) 
"""
Base.:+(bc::BoundConstraints,τ::AbstractArray) = BoundConstraints(lower_bound(bc).+τ,upper_bound(bc).+τ)
    
