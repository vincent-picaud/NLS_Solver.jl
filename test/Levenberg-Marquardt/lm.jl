@testset "LM algorithm" begin

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
    function eval_fobj(nls::Rosenbrock,θ::AbstractVector{<:Real})
        @assert length(θ)==parameter_size(nls)

        r1 = 1-θ[1]
        r2 = 10*(θ[2]-θ[1]^2)

        (r1*r1+r2*r2)/2
    end

    function eval_fobj_J(nls::Rosenbrock,θ::AbstractVector{<:Real})
        @assert length(θ)==parameter_size(nls)

        r1 = 1-θ[1]
        r2 = 10*(θ[2]-θ[1]^2)

        ∂1r1 = -1
        ∂2r1 =  0
        ∂1r2 =  -20*θ[1]
        ∂2r2 =  +10
        
        f = (r1*r1+r2*r2)/2
        J = typeof(f)[[ ∂1r1 ∂2r1 ] [ ∂1r2 ∂2r2 ]]
    end


    @testset "Rosenbrock def" begin
        fobj = Rosenbrock()
        
        @assert parameter_size(fobj) == 2

	θ=Float64[1,0]

        @test eval_fobj(fobj,θ) ≈ 50
        @test eval_fobj_J(fobj,θ) ≈ Float64[[-1 0] [-20 10]]
    end
    
end
 
