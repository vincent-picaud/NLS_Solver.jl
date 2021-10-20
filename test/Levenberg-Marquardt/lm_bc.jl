@testset "LM algorithm + bound constraints" begin

    @testset "Rosenbrock" begin
        using NLS_Solver: Levenberg_Marquardt
        
        nls = Rosenbrock()
        θ=Float64[-10;-10.5]
        bc=BoundConstraints(-40*ones(2),+20*ones(2))
        result=Levenberg_Marquardt_BC(nls, θ, bc)

        @test result
        # @test converged(result)
        # @test solution(result) ≈ Float64[1;1]
        # @test iteration_count(result) == 8
    end
        
end
 
