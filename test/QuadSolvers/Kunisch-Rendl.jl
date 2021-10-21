using LinearAlgebra: Symmetric, dot

@testset "Kunisch-Rendl" begin
    using Random: seed!
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


@testset "problem 1 uplow = $UL" for UL ∈ (:U,:L)
    using NLS_Solver: check_first_order
    
    n=5
    A=[Rational{Int}(1,i+j-1) for i in 1:n, j in 1:n]

    Q=Symmetric(A,UL)
    q=Rational{Int}[-1 for i in 1:n]
    bc=BoundConstraints(Int,n)
    x_init=zeros(Int,n)

    conf = Kunisch_Rendl_Conf()
    result = solve(Q,q,x_init,bc,conf)

    x_sol = solution(result)
    @test converged(result)
    @test iteration_count(result)==7
    @test objective_value(result) ≈ 1/2*dot(x_sol,Q*x_sol)+dot(q,x_sol)
    @test 1+check_first_order(Q,q,x_sol,bc) ≈ 1+0


    # test with a reg schedule
    conf_2 =  Kunisch_Rendl_Conf(reg_schedule=ExpRegularizationSchedule(factor=4,burning_last_iter=2))
    result_2 = solve(Q,q,x_init,bc,conf_2)

    x_sol_2 = solution(result_2)
    @test converged(result_2)
    @test iteration_count(result_2)==4 # note fewer iterations
    @test objective_value(result_2) ≈ 1/2*dot(x_sol_2,Q*x_sol_2)+dot(q,x_sol_2)
    @test 1+check_first_order(Q,q,x_sol_2,bc) ≈ 1+0

    # check that result = result_2
    @test x_sol ≈ x_sol_2
    @test objective_value(result) ≈ objective_value(result_2)
end

@testset "problem 2" begin
     using NLS_Solver: check_first_order
    
    Q=Symmetric(Float64[[30 20 15]
                        [20 15 12]
                        [15 12 10]],:L)

    q=-Float64[1:3;]
    bc=BoundConstraints(zeros(3),Float64[1:3;])
    x_init = zeros(3)

    conf = Kunisch_Rendl_Conf()
    result = solve(Q,q,x_init,bc,conf)

    x_sol = solution(result)
    τ = multiplier_τ(result)

    @test converged(result)
    @test iteration_count(result) == 4
    @test τ[1] ≈ -3.5
    @test τ[2] ≈ -1.6
    @test τ[3] == 0
    @test objective_value(result) ≈ 1/2*dot(x_sol,Q*x_sol)+dot(q,x_sol)
    @test 1+check_first_order(Q,q,x_sol,bc) ≈ 1+0

end 
