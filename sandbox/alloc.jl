# ################################################################
#
# CAVEAT: THIS EXAMPLE DOES NOT WORK YET
#
# REASON: the use of StaticArrays causes some problems:
#         https://github.com/JuliaArrays/StaticArrays.jl/issues/971
#        
# NOTE:   Rosenbrock_static is simply an alias of Rosenbrock and not a
#         static version of Rosenbrock as it should be.
#
# -> TO FIX
#
# ################################################################

# An example to count alloc static vs dynamic arrays
#
using Pkg, BenchmarkTools

Pkg.activate("../")

using NLS_Solver

include("../test/TestProblems/TestProblems.jl")

θ = [1.0,1.0]

obj = Rosenbrock()
obj_static = Rosenbrock_Static()

conf = LevenbergMarquardt_Conf()

println("Dynamic : ", @btime solve($obj, $θ, $conf))
println("Static  : ", @btime solve($obj_static, $θ, $conf))

# Sanity check
#
result = solve(obj, θ, conf)
result_static = solve(obj_static, θ, conf)

@assert converged(result)
@assert converged(result_static)
@assert solution(result) ≈ solution(result_static)

println("Sanity check ok")

# Constrained problem

bc = BoundConstraints(2)

conf = LevenbergMarquardt_BC_Conf()

println("Dynamic : ", @btime solve($obj, $θ, $bc, $conf))
println("Static  : ", @btime solve($obj_static, $θ, $bc, $conf))

#
# Debug runtime error for StaticArrays
#
# using Debugger
# breakpoint(joinpath(@__DIR__,"src/Levenberg-Marquardt/lm_bc.jl"),152)
# @enter solve(obj_static, θ, bc, conf)
