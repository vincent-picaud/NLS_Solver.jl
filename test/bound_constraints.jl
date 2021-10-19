@testset "BoundConstraints" begin

    @testset "Basic x F64, bc $T " for T in (Float64, Int)
        n=5
        bc = BoundConstraints(T,n)
        x = rand(n)
        
        @test lower_bound(bc)==zeros(T,n)
        @test upper_bound(bc)==ones(T,n)
        @test x ∈ bc
    end

    @testset "BC operations type=$T " for T in (Float64, Int)
        n=5
        bc = BoundConstraints(T,n)

        bc_copy = deepcopy(bc)
        
        bc = bc - ones(Int,n)
        
        @test lower_bound(bc)==-ones(T,n)
        @test upper_bound(bc)==zeros(T,n)
        @test eltype(bc)==T

        bc = bc + ones(T,n)

        @test lower_bound(bc)≈zeros(T,n)
        @test upper_bound(bc)≈ones(T,n)
    end
    
end 

