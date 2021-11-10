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

    grad=eval_nls_∇fobj(r,J)

    @test grad ≈ [200.0; -100.0]

    H=eval_nls_∇∇fobj(J)
    
    @test H ≈ Float64[401.0 -200.0; -200.0 100.0]
end

@testset "QuadModel wrapper" begin 

    using LinearAlgebra: Symmetric, dot

    nls = Rosenbrock()
    quad = NLSProblemAsQuadraticModel(nls)

    θ=Float64[1;0]

    (r, J) = eval_r_J(nls,θ)

    @test parameter_size(quad) == parameter_size(nls)
    @test eval_f(quad,θ) ≈ dot(r,r)/2
    @test all(eval_f_∇f(quad,θ) .≈ (dot(r,r)/2,J'*r))
    @test all(eval_f_∇f_∇∇f(quad,θ) .≈ (dot(r,r)/2,J'*r,Symmetric(J'*J)))
end
