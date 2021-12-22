@testset "LM algorithm + bound constraints" begin

    create_bc(l,u;n::Int) = BoundConstraints(l*ones(n),u*ones(n))
    
    @testset "Rosenbrock sol ∈ I" begin
        using NLS_Solver: LevenbergMarquardt_BC
        
        nls = Rosenbrock()
        θ=Float64[-10;-10]
        bc=create_bc(-1,2,n=parameter_size(nls))
        result=LevenbergMarquardt_BC(nls, θ, bc)

        @test converged(result)
        @test solution(result) ≈ Float64[1;1]
        @test iteration_count(result) == 14
    end

    @testset "Rosenbrock sol ∈ ∂I" begin
        using NLS_Solver: LevenbergMarquardt
        
        nls = Rosenbrock()
        θ=Float64[-10;-10]
        bc=create_bc(-1,1,n=parameter_size(nls))
        result=LevenbergMarquardt_BC(nls, θ, bc)

        @test converged(result)
        @test solution(result) ≈ Float64[1;1]
        @test iteration_count(result) == 14
    end

      @testset "Rosenbrock sol ∉ I" begin
        using NLS_Solver: LevenbergMarquardt
        
        nls = Rosenbrock()
        θ=Float64[-10;-10]
        bc=create_bc(-10,-2,n=parameter_size(nls))
        result=LevenbergMarquardt_BC(nls, θ, bc)

        @test converged(result)
        @test solution(result) ≈ Float64[-2;-2] # OK (MMA)
        @test iteration_count(result) == 4
    end

    # ================================================================
    
    @testset "Rosenbrock, use solve interface" begin
        nls = Rosenbrock()
        θ=Float64[-10;-10]
        bc=create_bc(-1,2,n=parameter_size(nls))

        conf = LevenbergMarquardt_BC_Conf()
        result=solve(nls, θ, bc, conf)

        @test converged(result)
        @test solution(result) ≈ Float64[1;1]
        @test iteration_count(result) == 14
    end

end
 
