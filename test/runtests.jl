using NLS_Solver
using Test

include("TestProblems/TestProblems.jl")

@testset "NLS_Solver.jl" begin

    @test length(detect_ambiguities(NLS_Solver))==0

    # ================================================================
    
    @testset "TestProblems" begin
        include("test_problems.jl")
    end

    # ================================================================
    
    include("bound_constraints.jl")
    include("regularization_schedule.jl")

    @testset "Quadsolvers" begin
        include("QuadSolvers/Kunisch-Rendl.jl")
    end

    @testset "Levenberg-Marquardt" begin
        include("Levenberg-Marquardt/lm.jl")
    end

end
