function single_run(example)
  alpha, beta, gamma, a, b, u, u_x, f = examples(example); ne = 2^3

  run_values = RunValues(alpha, beta, gamma, f, u)

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, ne)

  malha = monta_malha(ne, base, a, b)

  C, EQoLG, xPTne = solveSys(run_values, malha)

  return C
end

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

display(single_run(1) == correct_result_ex1)

display(single_run(2) == correct_result_ex2)



# P2 = [3, 7, 23]
# n_dim = length(P2)

# b2_1 = @benchmark PHI(P2, n_dim)
# b2_2 = @benchmark PHI_original(P2, n_dim)

# PHI(P2, n_dim)