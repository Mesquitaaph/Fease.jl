using MyProject, Test, BenchmarkTools

@testset "MyProject.jl" begin
    # Write your tests here.
    
    P1 = [3]
    n_dim = length(P1)
    
    @test PHI(P1, n_dim) == PHI_original(P1, n_dim)
    
    
    P2 = [3, 7]
    n_dim = length(P2)
    @test PHI(P2, n_dim) == PHI_original(P2, n_dim)

    # b2_1 = @benchmark PHI(P2, n_dim)
    # b2_2 = @benchmark PHI_original(P2, n_dim)

    # display(b2_1)
    # display(b2_2)
end
