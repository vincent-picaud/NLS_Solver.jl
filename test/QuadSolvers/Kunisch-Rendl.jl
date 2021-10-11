@testset "Kunisch-Rendl" begin
    using Random: seed!
    using LinearAlgebra: Symmetric
    using NLS_Solver: BoundContraintState_Enum, restrict_to_inactive!
    
    seed!(1234)

    n=10
    A=[Rational{Int}(1,i+j-1) for i in 1:n, j in 1:n]

    Z=rand(instances(BoundContraintState_Enum),n)
    lb=zeros(Rational{Int},n)
    ub=ones(Rational{Int},n)

    Q_result = Rational{Int64}[1//1 0//1 1//3 1//4 0//1 0//1 0//1 0//1 0//1 1//10; 0//1 1//1 0//1 0//1 0//1 0//1 0//1 0//1 0//1 0//1; 1//3 0//1 1//5 1//6 0//1 0//1 0//1 0//1 0//1 1//12; 1//4 0//1 1//6 1//7 0//1 0//1 0//1 0//1 0//1 1//13; 0//1 0//1 0//1 0//1 1//1 0//1 0//1 0//1 0//1 0//1; 0//1 0//1 0//1 0//1 0//1 1//1 0//1 0//1 0//1 0//1; 0//1 0//1 0//1 0//1 0//1 0//1 1//1 0//1 0//1 0//1; 0//1 0//1 0//1 0//1 0//1 0//1 0//1 1//1 0//1 0//1; 0//1 0//1 0//1 0//1 0//1 0//1 0//1 0//1 1//1 0//1; 1//10 0//1 1//12 1//13 0//1 0//1 0//1 0//1 0//1 1//19]
    q_result = Rational{Int64}[28//15, -1//1, 197//56, 1597//360, -1//1, -1//1, 0//1, 0//1, 0//1, 23629//2310]

    # Store upper part
    #
    Q=Symmetric(A,:U)
    q=Rational{Int}[i for i in 1:n]
    
    restrict_to_inactive!(Q, q, Z, lb, ub)

    @test Q==Q_result
    @test q==q_result

    # Store lower part
    #
    Q=Symmetric(A,:L)
    q=Rational{Int}[i for i in 1:n]

    restrict_to_inactive!(Q, q, Z, lb, ub)

    @test Q==Q_result
    @test q==q_result

end
