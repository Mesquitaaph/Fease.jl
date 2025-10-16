correct_result_ex1 = [0.018379439835948527
                      0.03508858369021478
                      0.048430969064432186
                      0.05665739050754408
                      0.057938496919442535
                      0.05033613330919126
                      0.031772978053653186]

correct_result_ex2 = [0.3831328931456573
                      0.7079372764182132
                      0.9249646268234876
                      1.0011744976201078
                      0.9249646268234876
                      0.7079372764182131
                      0.38313289314565724]

function single_run_1D(example)
  run_values = examples_1D(example)

  ne = 2^3
  a, b = 0, 1

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, ne)

  malha = monta_malha_1D_uniforme(base, ne, a, b)

  C = solve_sys_poisson(run_values, malha)

  # plot_solucao_aproximada(C, malha)

  display("Exemplo $example - Caso 1D")

  return C
end
display(single_run_1D(1) ≈ correct_result_ex1)
display(single_run_1D(2) ≈ correct_result_ex2)

function single_run_2D(ruido::Bool = false)
  run_values = examples_2D(1)

  Nx1, Nx2 = 4, 3
  ne = Nx1 * Nx2

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, ne)

  a = (0.0, 0.0)
  b = (1.0, 1.0)
  malha = monta_malha_2D_uniforme(base, Nx1, Nx2, a, b)

  correct_result_ex1 = [0.659197679603509
                        0.9322462987801567
                        0.659197679603509
                        0.6591976796035091
                        0.9322462987801565
                        0.659197679603509]

  correct_result_ex2 = [0.7065102277637031
                        0.969014464383587
                        0.6584488958892918
                        0.7346980208325475
                        0.8893466262809584
                        0.7082557097977525]

  if ruido
    # Exemplo 2: Adiciona ruído nos nós internos
    display("Exemplo 2: Adiciona ruído nos nós internos")
    Random.seed!(42) # Define uma semente para reprodutibilidade
    malha = malha2D_adiciona_ruido(malha)
    correct_result = correct_result_ex2
  else
    # Exemplo 1: Malha uniforme de retângulos
    display("Exemplo 1: Malha uniforme de retângulos")
    correct_result = correct_result_ex1
  end

  X₁, X₂ = malha.coords
  (; u) = run_values
  # Calcula a solução exata nos nós internos da malha
  C_exato = u.(X₁[2:(end-1), 2:(end-1)], X₂[2:(end-1), 2:(end-1)])
  C_exato = C_exato[:]

  K = monta_K_quadrilatero(run_values, malha)
  F = monta_F_quadrilatero(run_values, malha)
  C_2D = K \ F

  C = solve_sys_poisson(run_values, malha)

  display(C_2D ≈ correct_result)
  display(C ≈ correct_result)
  display(C_2D ≈ C)

  return C
end
single_run_2D(false)
single_run_2D(true)

function test_monta_F_1D()
  example = 1 # ou 2
  alpha, beta, gamma, a, b, u, u_x, f = examples_1D(example)
  ne = 2^3

  run_values = RunValues(alpha, beta, gamma, f, u)

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, ne)

  malha = monta_malha_1D_uniforme(base, ne, a, b)

  F_1D, xPTne = montaF_1D(run_values, malha)
  F_geral = montaF_geral(run_values, malha)

  display("Caso 1D: Exemplo 1")
  return display(F_1D ≈ F_geral)
end
# test_monta_F_1D()

function test_monta_F_2D()
  # Carrega os parâmetros de entrada da EDP
  alpha, beta, _, _, _, u, _, f = examples_2D(1)
  run_values = RunValues(alpha, beta, 0.0, f, u)

  # Define parâmetros da malha e monta a estrutura inicial
  Nx1, Nx2 = 4, 3

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, Nx1 * Nx2)

  a = (0.0, 0.0)
  b = (1.0, 1.0)
  malha = monta_malha_2D_uniforme(base, Nx1, Nx2, a, b)

  display("Exemplo 1: Malha uniforme de retângulos")
  F_2D = monta_F_quadrilatero(run_values, malha)
  F_geral = montaF_geral(run_values, malha)

  correct_F_2D = [0.9164924375447934
                  1.2961160349882237
                  0.9164924375447934
                  0.9164924375447936
                  1.296116034988224
                  0.9164924375447936]

  display(F_2D ≈ F_geral)

  # Exemplo 2: Malha com ruído nos nós internos
  display("Exemplo 2: Adiciona ruído nos nós internos")
  Random.seed!(42)  # Define uma semente para reprodutibilidade
  malha = malha2D_adiciona_ruido(malha)

  F_2D = monta_F_quadrilatero(run_values, malha)
  F_geral = montaF_geral(run_values, malha)

  # correct_F_2D = [
  #   0.9164924375447934
  #   1.2961160349882237
  #   0.9164924375447934
  #   0.9164924375447936
  #   1.296116034988224
  #   0.9164924375447936
  # ]

  return display(F_2D ≈ F_geral)
end
# test_monta_F_2D()

function test_monta_K_1D()
  example = 1
  alpha, beta, gamma, a, b, u, u_x, f = examples_1D(example)
  ne = 2^3

  run_values = RunValues(alpha, beta, gamma, f, u)

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, ne)

  malha = monta_malha_1D_uniforme(base, ne, a, b)

  K_1D = montaK_1D(run_values, malha)
  display(K_1D)
  K_geral = montaK_geral(run_values, malha)

  display("Caso 1D: Exemplo 1")
  return display(K_1D ≈ K_geral)
end
# test_monta_K_1D()

function test_monta_K_2D()
  # Carrega os parâmetros de entrada da EDP
  alpha, beta, _, _, _, u, _, f = examples_2D(1)
  run_values = RunValues(alpha, beta, 0.0, f, u)

  # Define parâmetros da malha e monta a estrutura inicial
  Nx1, Nx2 = 4, 3

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, Nx1 * Nx2)

  a = (0.0, 0.0)
  b = (1.0, 1.0)
  malha = monta_malha_2D_uniforme(base, Nx1, Nx2, a, b)

  # display("Exemplo 1: Malha uniforme de retângulos")
  # # display("monta_K_quadrilatero")
  # K_2D = monta_K_quadrilatero(run_values, malha)
  # # display(K_2D)
  # # display("montaK_geral")
  # K_geral = montaK_geral(run_values, malha)
  # # display(K_geral)

  # display(K_2D ≈ K_geral)

  # Exemplo 2: Malha com ruído nos nós internos
  display("Exemplo 2: Adiciona ruído nos nós internos")
  Random.seed!(42)  # Define uma semente para reprodutibilidade
  malha = malha2D_adiciona_ruido(malha)

  # display("monta_K_quadrilatero")
  K_2D = monta_K_quadrilatero(run_values, malha)
  # display(K_2D)
  # display("montaK_geral")
  K_geral = montaK_geral(run_values, malha)
  # display(K_geral)

  return display(K_2D ≈ K_geral)
end
# test_monta_K_2D()
