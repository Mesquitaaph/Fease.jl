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

function single_run_1D(example)
  alpha, beta, gamma, a, b, u, u_x, f = examples(example); ne = 2^3

  run_values = RunValues(alpha, beta, gamma, f, u)

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, ne)

  malha = monta_malha_1D_uniforme(ne, base, a, b)

  C, EQoLG, xPTne = solveSys_geral(run_values, malha)

  return C
end

function ≈(v1, v2)
  return all(1e-15 .> (v1 - v2))
end

# display(single_run_1D(1) ≈ correct_result_ex1)
# display(single_run_1D(2) ≈ correct_result_ex2)


function single_run_2D()
  # Carrega os parâmetros de entrada da EDP
  α, β, f, u = exemplo1()

  run_values = RunValues(α, β, 0.0, f, u)

	# Define parâmetros da malha e monta a estrutura inicial
  Nx1, Nx2 = 4, 3

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, Nx1*Nx2)
  
  a = (0.0, 0.0); b = (1.0, 1.0)
  malha = monta_malha_2D_uniforme(base, Nx1, Nx2, a ,b)

	display("Exemplo 1: Malha uniforme de retângulos")

  # Monta matriz K, vetor F e resolve o sistema linear Kc = F
  K = monta_K_quadrilatero(run_values, malha)
  F = monta_F_quadrilatero(run_values, malha)
  c = K \ F
  
  correct_result_ex1 = [
    0.659197679603509
    0.9322462987801567
    0.659197679603509
    0.6591976796035091
    0.9322462987801565
    0.659197679603509
  ]
  
  display(c ≈ correct_result_ex1)
  
	# Exemplo 2: Malha com ruído nos nós internos
	display("Exemplo 2: Adiciona ruído nos nós internos")
 	Random.seed!(42)  # Define uma semente para reprodutibilidade
	malha = malha2D_adiciona_ruido(malha)

  # Monta novamente K, F e resolve o sistema com a malha perturbada
  K = monta_K_quadrilatero(run_values, malha)
  F = monta_F_quadrilatero(run_values, malha)
  c = K \ F

  correct_result_ex2 = [
    0.7065102277637031
    0.969014464383587
    0.6584488958892918
    0.7346980208325475
    0.8893466262809584
    0.7082557097977525
  ]
    
  display(c ≈ correct_result_ex2)
end

single_run_2D()