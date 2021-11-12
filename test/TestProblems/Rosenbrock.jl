
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

function eval_r(nls::Rosenbrock,θ::AbstractVector{T}) where T
    @assert length(θ)==parameter_size(nls)

    r = T[1-θ[1],
          10*(θ[2]-θ[1]^2)]
end

function eval_r_J(nls::Rosenbrock,θ::AbstractVector{T}) where T
    @assert length(θ)==parameter_size(nls)

    r = SVector{2,T}(1-θ[1],
                     10*(θ[2]-θ[1]^2))

    # CAVEAT: SMatrix are filled column by column
    #
    # The Jacobian
    #
    # | ∂1r1, ∂2r1 |
    # | ∂1r2, ∂2r2 |
    #
    # must be created by
    #
    # SMatrix(∂1r1,∂1r2, # col 1
    #         ∂2r1,∂2r2) # col 2
    #
    # which is not intuitive and error prone.
    #
    # The alternative to is to use the @SMatrix macro
    # that follows the "natural" order (row by row)
    #
    J = T[      -1   +0;     # ∂1r1, ∂2r1
          -20*θ[1]  +10]     # ∂1r2, ∂2r2

    (r, J)
end

# ================================================================

# Define the same Rosenbrock function but using Statics Arrays... This
# cause some problems: https://github.com/JuliaArrays/StaticArrays.jl/issues/971
#
struct Rosenbrock_Static <: AbstractNLS
end


parameter_size(::Rosenbrock_Static) = 2
residue_size(::Rosenbrock_Static) = 2

function eval_r(nls::Rosenbrock_Static,θ::AbstractVector{T}) where T
    @assert length(θ)==parameter_size(nls)

    r = T[1-θ[1],
          10*(θ[2]-θ[1]^2)]
end


function eval_r(nls::Rosenbrock_Static,θ::AbstractVector{T}) where T
    @assert length(θ)==parameter_size(nls)

    r = @SVector T[1-θ[1],
                   10*(θ[2]-θ[1]^2)]
end

function eval_r_J(nls::Rosenbrock_Static,θ::AbstractVector{T}) where T
    @assert length(θ)==parameter_size(nls)

    r = SVector{2,T}(1-θ[1],
                     10*(θ[2]-θ[1]^2))

    # CAVEAT: SMatrix are filled column by column
    #
    # The Jacobian
    #
    # | ∂1r1, ∂2r1 |
    # | ∂1r2, ∂2r2 |
    #
    # must be created by
    #
    # SMatrix(∂1r1,∂1r2, # col 1
    #         ∂2r1,∂2r2) # col 2
    #
    # which is not intuitive and error prone.
    #
    # The alternative to is to use the @SMatrix macro
    # that follows the "natural" order (row by row)
    #
    J = @SMatrix T[      -1   +0;     # ∂1r1, ∂2r1
                    -20*θ[1]  +10]     # ∂1r2, ∂2r2

    (r, J)
end
