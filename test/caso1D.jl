function single_run_1D(example)
  run_values = examples_1D(example)
  ne = 2^3

  a, b = 0, 1

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, ne)

  malha = monta_malha_1D_uniforme(base, ne, a, b)

  C = solve_sys_poisson(run_values, malha)

  return C
end

@testset "caso1D.jl" begin
  correct_result_ex1 = [0.018379439835948527
                        0.03508858369021478
                        0.048430969064432186
                        0.05665739050754408
                        0.057938496919442535
                        0.05033613330919126
                        0.031772978053653186]
  @test single_run_1D(1) ≈ correct_result_ex1

  correct_result_ex2 = [0.3831328931456573
                        0.7079372764182132
                        0.9249646268234876
                        1.0011744976201078
                        0.9249646268234876
                        0.7079372764182131
                        0.38313289314565724]
  @test single_run_1D(2) ≈ correct_result_ex2
end
