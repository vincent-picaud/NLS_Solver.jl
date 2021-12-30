export Abstract_BC_Solver_Result

# Result abstraction for bound constrained problems ================
#
# For the moment only use an alias as we do not have specific methods
#
@doc raw"""

The structure returned by [`solve`](@ref) when using the
[`LevenbergMarquardt_BC_Conf`](@ref) method.

See [`Abstract_Solver_Result`](@ref) 
"""
const Abstract_BC_Solver_Result = Abstract_Solver_Result

