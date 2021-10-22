@doc raw"""

A structure to store [`Levenberg_Marquardt_BC`](@ref) results. 

This structure subtypes [`LevenbergMarquardt_BC_Result`](@ref)

"""
const LevenbergMarquardt_BC_Result = LevenbergMarquardt_Result{T} where {T<:Real}
