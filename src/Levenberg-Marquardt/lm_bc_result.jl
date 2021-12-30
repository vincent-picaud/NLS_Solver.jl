export LevenbergMarquardt_BC_Result 

@doc raw"""

```julia
struct LevenbergMarquardt_BC_Result{T<:Real} <:  Abstract_BC_Solver_Result
    ...
end 
```

This structure subtypes [`Abstract_BC_Solver_Result`](@ref)

"""
const LevenbergMarquardt_BC_Result = LevenbergMarquardt_Result{T} where {T<:Real}
