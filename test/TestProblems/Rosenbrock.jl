
# Define Rosenbrock function: https://en.wikipedia.org/wiki/Rosenbrock_function
#
# This function can be interpreted as a NLS pb:
# r_1 = (1-θ₁)
# r_2 = 10(θ₂-θ₁²)
#
# Jacobian is: [[-1, 0],[-20, 10]]
#
# Minimum for θ=(1,1) where f=1/2 r^2 = 0
#
# TODO: implement even n (here 2) generalization of the wiki page.
#
struct Rosenbrock <: AbstractNLS
end


parameter_size(::Rosenbrock) = 2
residue_size(::Rosenbrock) = 2

function eval_r!(r::AbstractVector,
                 nls::Rosenbrock,θ::AbstractVector)
    @assert length(θ)==parameter_size(nls)
    @assert length(r)==residue_size(nls)

    r[1] = 1-θ[1]
    r[2] = 10*(θ[2]-θ[1]^2)

    r
end

function eval_r_J!(r::AbstractVector,J::AbstractMatrix,
                   nls::Rosenbrock,θ::AbstractVector)
    @assert length(θ)==parameter_size(nls)
    @assert length(r)==residue_size(nls)
    @assert size(J)==(parameter_size(nls),residue_size(nls))

    r[1] = 1-θ[1]
    r[2] = 10*(θ[2]-θ[1]^2)

    J[1,1] = -1        # ∂1r1
    J[1,2] =  0        # ∂2r1
    J[2,1] =  -20*θ[1] # ∂1r2
    J[2,2] =  +10      # ∂2r2

    (r, J)
end
