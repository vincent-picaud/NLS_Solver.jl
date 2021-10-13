@testset "BoundConstraints" begin
    n=5
    bc = BoundConstraints(n)
    x = rand(n)

    @test lower_bound(bc)==zeros(n)
    @test lower_bound(bc)==ones(n)
    @test x ∈ bc
end 
