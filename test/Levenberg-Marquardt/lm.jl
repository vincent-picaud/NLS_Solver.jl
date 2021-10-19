@testset "LM algorithm" begin

    @testset "LM Rosenbrock" begin
        nls = Rosenbrock()
        θ=Float64[1;0]

        result=Levenberg_Marquardt(nls, θ)

        @test converged(result)
        @test solution(result) ≈ Float64[1;1]
        @test iteration_count(result) == 8
    end
        
end
 
