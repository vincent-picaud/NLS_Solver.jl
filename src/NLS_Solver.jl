module NLS_Solver

include("regularization_schedule.jl")
include("bound_constraints.jl")

include("QuadSolvers/QuadSolvers.jl")

include("abstract_nls.jl")
include("abstract_nls_result.jl")
include("abstract_nls_conf.jl")
include("abstract_nls_solve.jl")

include("Levenberg-Marquardt/Levenberg-Marquardt.jl")

end
