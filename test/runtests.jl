using NLS_Solver
using Test

@testset "NLS_Solver.jl" begin

    @testset "Quadsolvers" begin
        include("QuadSolvers/Kunisch-Rendl.jl")
    end

end
