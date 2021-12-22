@testset "LM algorithm" begin

    @testset "Rosenbrock" for TYPE in (Rosenbrock, Rosenbrock_Static)
        using NLS_Solver: Levenberg_Marquardt
        
        nls = TYPE()
        θ=Float64[1;0]

        result=Levenberg_Marquardt(nls, θ)

        @test converged(result)
        @test solution(result) ≈ Float64[1;1]
        @test iteration_count(result) == 8
    end
    
    @testset "Rosenbrock with solve()" for TYPE in (Rosenbrock, Rosenbrock_Static)
        nls = TYPE()
        conf = LevenbergMarquardt_Conf()
        θ=Float64[1;0]

        result=solve( nls, θ, conf)

        @test converged(result)
        @test solution(result) ≈ Float64[1;1]
        @test iteration_count(result) == 8
    end

    @testset "Rosenbrock with solve() + distant initial point" for TYPE in (Rosenbrock, Rosenbrock_Static) 
        nls = TYPE()
        conf = LevenbergMarquardt_Conf()
        θ=Float64[-10;-10.5]

        result=solve( nls, θ, conf)

        @test converged(result)
        @test solution(result) ≈ Float64[1;1]
        @test iteration_count(result) == 30
    end

    @testset "PowellSingular" begin 
        nls = PowellSingular()
        conf = LevenbergMarquardt_Conf(ε_grad_inf_norm=1e-15,
                                        ε_step_2_norm=1e-15)
        θ=Float64[3;-1;0;1]
        result=solve( nls, θ, conf)

        @test converged(result)
        @test isapprox(solution(result),zeros(4),atol=1e-4)
        @test iteration_count(result) == 23
    end 
end
 
