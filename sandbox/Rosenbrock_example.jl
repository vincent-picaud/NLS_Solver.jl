using Revise
using NLS_Solver

struct Rosenbrock <: NLS_Solver.AbstractNLS
end

import NLS_Solver: parameter_size, residue_size, eval_r, eval_r_J

NLS_Solver.parameter_size(::Rosenbrock) = 2
NLS_Solver.residue_size(::Rosenbrock) = 2

function NLS_Solver.eval_r(nls::Rosenbrock,θ::AbstractVector{T}) where T
    @assert length(θ)==parameter_size(nls)

    sqrt(2)* T[ 1-θ[1], 10*(θ[2]-θ[1]^2) ]
end

function NLS_Solver.eval_r_J(nls::Rosenbrock,θ::AbstractVector{T}) where T
    @assert length(θ)==parameter_size(nls)

    r = sqrt(2)* T[ 1-θ[1], 10*(θ[2]-θ[1]^2) ]
    J = sqrt(2)* T[ -1 0; -20*θ[1] 10]

    (r,J)
end


nls = create_NLS_problem_using_ForwardDiff(2 => 2) do θ
  sqrt(2)* [ 1-θ[1], 10*(θ[2]-θ[1]^2) ]
end


conf = Levenberg_Marquardt_Conf()

θ_init = zeros(2)

result = solve(nls, θ_init, conf)
