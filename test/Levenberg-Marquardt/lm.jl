@testset "LM algorithm" begin


    @testset "Rosenbrock def" begin
        nls = Rosenbrock()
        
        @assert parameter_size(nls) == 2

	θ=Float64[1,0]

        (r, J) = eval_r_J(nls,θ)
        
        @test r ≈ Float64[0; -10]
        @test J ≈ Float64[[-1 0]; [-20 10]]
        @test eval_nls_fobj(r) ≈ 50

        grad=similar(r)

        eval_nls_∇fobj!(grad,r,J)
        @test grad ≈ [0; -100.0]
    end
    
end
 
