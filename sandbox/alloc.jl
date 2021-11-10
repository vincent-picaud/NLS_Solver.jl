# An example to count alloc static vs dynamic arrays
#
using Pkg, BenchmarkTools

Pkg.activate("../")

using NLS_Solver

include("../test/TestProblems/TestProblems.jl")

θ = [1.0,1.0]

obj = Rosenbrock()
obj_static = Rosenbrock_Static()

conf = Levenberg_Marquardt_Conf()

println("Dynamic : ", @ballocated solve($obj, $θ, $conf))
println("Static  : ", @ballocated solve($obj_static, $θ, $conf))

# Sanity check
#
result = solve(obj, θ, conf)
result_static = solve(obj_static, θ, conf)

@assert converged(result)
@assert converged(result_static)
@assert solution(result) ≈ solution(result_static)

println("Sanity check ok")
