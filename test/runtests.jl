using NLS_Solver
using Test

@testset "NLS_Solver.jl" begin

    include("bound_constraints.jl")
    include("regularization_schedule.jl")

    @testset "Quadsolvers" begin
        include("QuadSolvers/Kunisch-Rendl.jl")
    end

    @testset "Levenberg-Marquardt" begin
        include("Levenberg-Marquardt/lm.jl")
    end

end
