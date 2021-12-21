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

    @testset "Quadsolvers Directory" begin
        include("QuadSolvers/misc.jl")
        include("QuadSolvers/Kunisch-Rendl.jl")
    end

    @testset "Levenberg-Marquardt" begin
        include("Levenberg-Marquardt/lm.jl")
    end

    @testset "Bound constrained Levenberg-Marquardt" begin
        include("Levenberg-Marquardt/lm_bc.jl")
    end

end
