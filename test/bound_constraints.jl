@testset "BoundConstraints" begin
    n=5
    bc = BoundConstraints(n)
    x = rand(n)

    @test lower_bound(bc)==zeros(Int,n)
    @test lower_bound(bc)==ones(Int,n)
    @test x ∈ bc
end 
