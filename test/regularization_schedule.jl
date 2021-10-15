@testset "RegularizationSchedule" begin

    @testset "No" begin
        rs = NoRegularizationSchedule()

        @test burning_phase(rs,1) == false
        @test burning_phase(rs,4) == false
        @test burning_phase(rs,5) == false
        @test burning_phase(rs,6) == false

        @test regularization_factor(rs,1) ≈ 1
        @test regularization_factor(rs,4) ≈ 1
        @test regularization_factor(rs,5) ≈ 1
        @test regularization_factor(rs,6) ≈ 1
    end
    
    @testset "Exp" begin
        rs = ExpRegularizationSchedule(factor=10,burning_last_iter=5)

        @test burning_phase(rs,1) == true
        @test burning_phase(rs,4) == true
        @test burning_phase(rs,5) == true
        @test burning_phase(rs,6) == false
        @test burning_phase(rs,7) == false

        @test regularization_factor(rs,1) ≈ 10
        @test regularization_factor(rs,4) > 1
        @test regularization_factor(rs,5) > 1
        @test regularization_factor(rs,6) ≈ 1
        @test regularization_factor(rs,7) ≈ 1
    end
    
end 

