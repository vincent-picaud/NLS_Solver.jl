@testset "BoundConstraints" begin
    bc = BoundConstraints(5)
    x = rand(5)

    @test x ∈ bc
end 
