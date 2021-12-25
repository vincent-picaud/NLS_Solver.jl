export solve

using LinearAlgebra: Symmetric

#
# Generic interface to solve a quadratic optimization problem of the form:
#
#  min  1/2 <Qx,x> + <q,x>
# l≤x≤u
#
solve(Q::Symmetric{<:Real},
      q::AbstractVector{<:Real},
      x_init::AbstractVector{<:Real},
      bc::BoundConstraints{<:Real,1},
      conf::Abstract_BC_QuadSolver_Conf) = @assert(false,"To implement")

               
