@testset "LM algorithm" begin


    @testset "Rosenbrock def" begin
        fobj = Rosenbrock()
        
        @assert parameter_size(fobj) == 2

	θ=Float64[1,0]

        @test eval_fobj(fobj,θ) ≈ 50
        @test eval_fobj_J(fobj,θ) ≈ Float64[[-1 0] [-20 10]]
    end
    
end
 
