@testset "LM algorithm" begin

    @testset "LM Rosenbrock" begin
        using NLS_Solver: Levenberg_Marquardt
        
        nls = Rosenbrock()
        θ=Float64[1;0]

        result=Levenberg_Marquardt(nls, θ)

        @test converged(result)
        @test solution(result) ≈ Float64[1;1]
        @test iteration_count(result) == 8
    end
    
    @testset "LM Rosenbrock with solve()" begin
        
        nls = Rosenbrock()
        conf = Levenberg_Marquardt_Conf()
        θ=Float64[1;0]

        result=solve(conf, nls, θ)

        @test converged(result)
        @test solution(result) ≈ Float64[1;1]
        @test iteration_count(result) == 8
    end
        
end
 
