# ================================================================
# Damping factor
# ================================================================
#
@doc raw"""
```julia
Abstract damping factor
```

The interface is:
- [`get_damping_factor(::AbstractDampingFactor)`](@ref)
- [`increase_damping_factor(::AbstractDampingFactor)`](@ref)
- [`decrease_damping_factor(::AbstractDampingFactor)`](@ref)
- [`reset(::AbstractDampingFactor, damping_factor::Real)`](@ref)

Concrete impementations are:
- [`DampingFactor`](@ref)
"""
abstract type AbstractDampingFactor end

@doc raw"""
```julia
get_damping_factor(::AbstractDampingFactor)
```
Return the damping factor (often noted as μ)
"""
get_damping_factor(::AbstractDampingFactor) = @assert(false,"To implement")

@doc raw"""
```julia
increase_damping_factor(::AbstractDampingFactor)
```
Increase the damping factor (often noted as μ)
"""
increase_damping_factor(::AbstractDampingFactor) = @assert(false,"To implement")

@doc raw"""
```julia
decrease_damping_factor(::AbstractDampingFactor)
```
Decrease the damping factor (often noted as μ)
"""
decrease_damping_factor(::AbstractDampingFactor) = @assert(false,"To implement")

@doc raw"""
```julia
reset(::AbstractDampingFactor, damping_factor::Real)
```
Reset the damping factor (often noted as μ)
"""
reset(::AbstractDampingFactor, damping_factor::Real) = @assert(false,"To implement")

@doc raw"""
```julia
DampingFactor(damping_factor::Real)
```
Concrete implementation of [`AbstractDampingFactor`](@ref)

"""
struct DampingFactor <: AbstractDampingFactor
    _damping_factor::Float64
    _growing_factor::Float64
    _growing_factor_factor::Float64

    function DampingFactor(damping_factor::Real,
                           growing_factor::Real,
                           growing_factor_factor::Real)
        @assert damping_factor ≥ 0
        @assert growing_factor > 1
        @assert growing_factor_factor > 1

        new(damping_factor,growing_factor,growing_factor_factor)
    end
    
    function DampingFactor(damping_factor::Real)
        new(damping_factor,2,2)
    end 
                  
end 

get_damping_factor(dp::DampingFactor) = dp._damping_factor
function increase_damping_factor(dp::DampingFactor)
    DampingFactor(dp._damping_factor*dp._growing_factor,
                  dp._growing_factor*dp._growing_factor_factor,
                  dp._growing_factor_factor)
end
function decrease_damping_factor(dp::DampingFactor)
    growing = max(dp._growing_factor_factor,dp._growing_factor/dp._growing_factor_factor)
    DampingFactor(dp._damping_factor/growing,
                  growing,
                  dp._growing_factor_factor)
end

reset(dp::DampingFactor, damping_factor::Real) = DampingFactor(damping_factor,
                                                               dp._growing_factor_factor,
                                                               dp._growing_factor_factor)


# ****************************************************************

@doc raw"""
```julia
Abstract damping factor
```

The interface is:
- [`get_damping_factor(::AbstractDynamicDampingFactor)`](@ref)
- [`update_damping_factor(::AbstractDynamicDampingFactor,gain_ratio::Real)`](@ref)
- [`reset(::AbstractDynamicDampingFactor, damping_factor::Real)`](@ref)

Concrete impementations are:
- [`DynamicDampingFactor`](@ref)
"""
abstract type AbstractDynamicDampingFactor end

@doc raw"""
```julia
get_damping_factor(::AbstractDynamicDampingFactor)
```
Return the damping factor (often noted as μ)
"""
get_damping_factor(::AbstractDynamicDampingFactor) = @assert(false,"To implement")

@doc raw"""
```julia
update_damping_factor(::AbstractDynamicDampingFactor,gain_ratio::Real)
```
Update the damping factor (often noted as μ) according to the gain ration (often noted as ρ)
"""
update_damping_factor(::AbstractDynamicDampingFactor,gain_ratio::Real) = @assert(false,"To implement")

@doc raw"""
```julia
reset(::AbstractDynamicDampingFactor, damping_factor::Real)
```
Reset the damping factor (often noted as μ)
"""
reset(::AbstractDynamicDampingFactor,damping_μ::Real) = @assert(false,"To implement")

@doc raw"""
```julia
DynamicDampingFactor(damping_factor::Real)
```
Concrete implementation of [`AbstractDynamicDampingFactor`](@ref)

"""
struct DynamicDampingFactor <: AbstractDynamicDampingFactor
    _dp::DampingFactor

    function DynamicDampingFactor(damping_factor::Real)
        new(DampingFactor(damping_factor))
    end 
    function DynamicDampingFactor(dp::DampingFactor)
        new(dp)
    end 
                  
end 

get_damping_factor(dp::DynamicDampingFactor) = get_damping_factor(dp._dp)
function update_damping_factor(dp::DynamicDampingFactor,gain_ratio::Real)
    if gain_ratio>0
        # decrease damping_factor
        #
        μ = get_damping_factor(dp) * max(1/3,1-(2*gain_ratio-1)^3)
        return reset(dp,μ)
    end

    DynamicDampingFactor(increase_damping_factor(dp._dp))
end 
function reset(dp::DynamicDampingFactor,damping_factor::Real)
    DynamicDampingFactor(damping_factor)
end 
