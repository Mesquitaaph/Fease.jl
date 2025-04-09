using MyProject, Test, BenchmarkTools

@testset "MyProject.jl" begin
    # Write your tests here.
    
    # P1 = [3]
    # n_dim = length(P1)
    
    # @test PHI(P1, n_dim) == PHI_original(P1, n_dim)
    
    
    # P2 = [3, 7]
    # n_dim = length(P2)
    # @test PHI(P2, n_dim) == PHI_original(P2, n_dim)

    # @test PHI_original(P2, 3) == Nothing

    correct_result_ex1 = [
        0.018379439835948527
        0.03508858369021478
        0.048430969064432186
        0.05665739050754408
        0.057938496919442535
        0.05033613330919126
        0.031772978053653186
    ]

    correct_result_ex2 = [
        0.3831328931456573
        0.7079372764182132
        0.9249646268234876
        1.0011744976201078
        0.9249646268234876
        0.7079372764182131
        0.38313289314565724
    ]

    @test single_run_1D(1) ≈ correct_result_ex1

    @test single_run_1D(2) ≈ correct_result_ex2
end