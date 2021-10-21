@testset "LM algorithm + bound constraints" begin

    create_bc(l,u;n::Int) = BoundConstraints(l*ones(n),u*ones(n))
    
    @testset "Rosenbrock sol ∈ I" begin
        using NLS_Solver: Levenberg_Marquardt
        
        nls = Rosenbrock()
        θ=Float64[-10;-10]
        bc=create_bc(-1,2,n=parameter_size(nls))
        result=Levenberg_Marquardt_BC(nls, θ, bc)

        @test converged(result)
        @test solution(result) ≈ Float64[1;1]
        @test iteration_count(result) == 10
    end

    @testset "Rosenbrock sol ∈ ∂I" begin
        using NLS_Solver: Levenberg_Marquardt
        
        nls = Rosenbrock()
        θ=Float64[-10;-10]
        bc=create_bc(-1,1,n=parameter_size(nls))
        result=Levenberg_Marquardt_BC(nls, θ, bc)

        @test converged(result)
        @test solution(result) ≈ Float64[1;1]
        @test iteration_count(result) == 10
    end

      @testset "Rosenbrock sol ∉ I" begin
        using NLS_Solver: Levenberg_Marquardt
        
        nls = Rosenbrock()
        θ=Float64[-10;-10]
        bc=create_bc(-10,-2,n=parameter_size(nls))
        result=Levenberg_Marquardt_BC(nls, θ, bc)

        @test converged(result)
        @test solution(result) ≈ Float64[-2;-2] # OK (MMA)
        @test iteration_count(result) == 4
    end
        
end
 
