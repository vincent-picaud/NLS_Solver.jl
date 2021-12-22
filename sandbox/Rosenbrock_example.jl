using Revise

# The most straightforward approach ****************
#
using NLS_Solver

nls = create_NLS_problem_using_ForwardDiff(2 => 2) do θ
  sqrt(2)* [ 1-θ[1], 10*(θ[2]-θ[1]^2) ]
end


conf = LevenbergMarquardt_Conf()

θ_init = zeros(2)

result = solve(nls, θ_init, conf)

@assert converged(result)

θ_solution = solution(result)


# Same example without the do...end syntax ****************
#
Rosenbrock(θ::AbstractVector{T}) where T = sqrt(2)* T[ 1-θ[1], 10*(θ[2]-θ[1]^2) ]

nls = create_NLS_problem_using_ForwardDiff(θ->Rosenbrock(θ),2 => 2);
 
result = solve(nls, θ_init, conf)

@assert converged(result)

θ_solution = solution(result)

# Sub-typing AbstractNLS in details ****************
#
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

nls = Rosenbrock()

result = solve(nls, θ_init, conf)

@assert converged(result)

θ_solution = solution(result)

# Bound constrained problems ****************
#
conf = LevenbergMarquardt_BC_Conf()

θl = Float64[2,2]
θu = Float64[4,4]

bc = BoundConstraints(θl,θu)

result = solve(nls, θ_init, bc, conf)

@assert converged(result)

θ_solution = solution(result)


#
# Rosenbrock(θ::AbstractVector{T}) where T = sqrt(2)* T[ 1-θ[1], 10*(θ[2]-θ[1]^2) ]
#
# nls = create_NLS_problem_using_ForwardDiff(θ->Rosenbrock(θ),2 => 2);
 
