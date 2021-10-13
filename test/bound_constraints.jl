@testset "BoundConstraints" begin

    @testset "Basic x F64, bc $T " for T in (Float64, Int)
        n=5
        bc = BoundConstraints(T,n)
        x = rand(n)
        
        @test lower_bound(bc)==zeros(T,n)
        @test upper_bound(bc)==ones(T,n)
        @test x ∈ bc
    end
    
end 

