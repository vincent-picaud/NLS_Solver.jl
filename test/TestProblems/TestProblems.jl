# CAVEAT: do not forget this!
#
# -> TestProblem does not belong to the NLS_Solver module.
#    But one must overload these functions, hence we import them.
#
import NLS_Solver: parameter_size, residue_size, eval_r!, eval_r_J!

include("Rosenbrock.jl")
include("PowellSingular.jl")
