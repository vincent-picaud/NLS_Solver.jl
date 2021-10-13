@testset "Kunisch-Rendl" begin
    using Random: seed!
    using LinearAlgebra: Symmetric
    using NLS_Solver: BoundConstraint_Enum, restrict_to_inactive!
    
    seed!(1234)

    n=10
    A=[Rational{Int}(1,i+j-1) for i in 1:n, j in 1:n]

    Z=rand(instances(BoundConstraint_Enum),n)
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


@testset "problem 1" begin
    using NLS_Solver: Kunisch_Rendl,check_first_order,create_damping_schedule_nothing
    
    n=5
    A=[Rational{Int}(1,i+j-1) for i in 1:n, j in 1:n]

    Q=Symmetric(A,:U)
    q=Rational{Int}[-1 for i in 1:n]
    bc=BoundConstraints(Int,n)
    x_init=zeros(Int,n)

    (cv,x_sol,τ_sol) = Kunisch_Rendl(Q,q,x_init,bc,10,create_damping_schedule_nothing())

    @test 1+check_first_order(Q,q,x_sol,bc) ≈ 1+0
end

@testset "problem 2" begin
    using NLS_Solver: Kunisch_Rendl,check_first_order,create_damping_schedule_nothing

    Q=Symmetric(Float64[[30 20 15]
                        [20 15 12]
                        [15 12 10]],:L)

    q=-Float64[1:3;]
    bc=BoundConstraints(zeros(3),Float64[1:3;])
    x_init = zeros(3)
    (cv,x_sol,τ_sol) = Kunisch_Rendl(Q,q,x_init,bc,10,create_damping_schedule_nothing())

    @test 1+check_first_order(Q,q,x_sol,bc) ≈ 1+0
end 
