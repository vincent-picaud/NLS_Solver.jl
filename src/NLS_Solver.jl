module NLS_Solver
include("log_format.jl")
include("regularization_schedule.jl")
include("bound_constraints.jl")

include("QuadSolvers/QuadSolvers.jl")

include("abstract_nls.jl")
include("abstract_nls_forwarddiff.jl")

include("abstract_solver_result.jl")
include("abstract_bc_solver_result.jl")

include("abstract_solver_conf.jl")
include("abstract_bc_solver_conf.jl")

include("solver_interface.jl")
include("bc_solver_interface.jl")

include("Levenberg-Marquardt/Levenberg-Marquardt.jl")

end
