@testset "LM algorithm" begin
    using LinearAlgebra: Symmetric
    @testset "Rosenbrock def" begin
        nls = Rosenbrock()
        
        @assert parameter_size(nls) == 2
        @assert residue_size(nls) == 2

	θ=Float64[1;0]

        (r, J) = eval_r_J(nls,θ)
        
        @test r ≈ Float64[0; -10]
        @test J ≈ Float64[-1 0; -20 10]
        @test eval_nls_fobj(r) ≈ 50

        grad=similar(r)

        eval_nls_∇fobj!(grad,r,J)
        @test grad ≈ [200.0; -100.0]

        H=Symmetric(Matrix{Float64}(undef,parameter_size(nls),parameter_size(nls)))
        eval_nls_∇∇fobj!(H,J)
        @test H ≈ Float64[401.0 -200.0; -200.0 100.0]
    end

    @testset "LM Rosenbrock" begin
        nls = Rosenbrock()
        θ=Float64[1;0]

        status=Levenberg_Marquardt(nls, θ)

        @test status
    end
        
end
 
