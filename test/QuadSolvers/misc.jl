
@testset "Q+μD" begin
    using LinearAlgebra: Symmetric, I
    using NLS_Solver: compute_Q_μD
    
    n = 5
    
    ref = 2.0*I(n)

    Q=Symmetric(zeros(n,n))

    Q1 = compute_Q_μD(Q,2.0,I)
    @test Q1 ≈ ref

    Q2 = compute_Q_μD(Q,2.0,ones(n))
    @test Q2 ≈ ref
end
 
