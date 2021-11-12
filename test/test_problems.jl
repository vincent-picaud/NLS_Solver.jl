@testset "Rosenbrock" begin
    using LinearAlgebra: Symmetric
    
    nls = Rosenbrock()
    
    @test parameter_size(nls) == 2
    @test residue_size(nls) == 2

    θ=Float64[1;0]

    (r, J) = eval_r_J(nls,θ)
    
    @test r ≈ Float64[0; -10]
    @test J ≈ Float64[-1 0; -20 10]
    @test eval_nls_fobj(r) ≈ 50

    grad = eval_nls_∇fobj(r,J)
    @test grad ≈ [200.0; -100.0]

    H = eval_nls_∇∇fobj(J)
    @test H ≈ Float64[401.0 -200.0; -200.0 100.0]
end
