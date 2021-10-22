module NLS_Solver
include("log_format.jl")
include("regularization_schedule.jl")
include("bound_constraints.jl")

include("QuadSolvers/QuadSolvers.jl")

include("abstract_nls.jl")

include("abstract_nls_result.jl")
include("abstract_nls_bc_result.jl")

include("abstract_nls_conf.jl")
include("abstract_nls_conf_bc.jl")

include("abstract_nls_solve.jl")
include("abstract_nls_solve_bc.jl")

include("Levenberg-Marquardt/Levenberg-Marquardt.jl")

end
