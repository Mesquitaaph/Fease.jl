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

function single_run(example)
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

display(single_run(1) ≈ correct_result_ex1)
display(single_run(2) ≈ correct_result_ex2)

function test_montaLG()
  ne = 5

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, ne)
  LG_1D = montaLG(ne, base)
  println("Teste montagem LG\n")

  println("LG 1D")
  display(LG_1D == montaLG_geral(5))

  println("LG 2D")
  display(monta_LG_2D(5, 4) == montaLG_geral(5, 4))
end
# test_montaLG()

function test_montaEQ()
  ne = 5

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, ne)
  neq = base.neq
  EQ_1D = montaEQ(ne, neq, base)
  println("\nTeste montagem EQ\n")

  println("EQ 1D")  
  display(EQ_1D == montaEQ_geral(5))
  # display(EQ_1D); display(montaEQ_geral(5))

  println("\nEQ 2D")
  display(monta_EQ_2D(5, 4)[2] == montaEQ_geral(5, 4))
  # display(monta_EQ_2D(5, 4)); display(montaEQ_geral(5, 4))
  for nx1 in 1:50
    for nx2 in 1:50
      if !(monta_EQ_2D(nx1, nx2)[2] == montaEQ_geral(nx1, nx2))
        display("Falha")
      end
    end
  end
end
# test_montaEQ()

function test_ϕ()
  println("Teste da ϕ")
  P, W = legendre(10)

  print("ϕ 1D: ")
  teste_1D = true

  for ξ₁ in P
    if !(ϕ_1D(ξ₁) == ϕ_geral(ξ₁))
      teste_1D = false
      break
    end
  end
  println(teste_1D ? "ok" : "nok")

  print("ϕ 2D: ")
  teste_2D = true

  for ξ₁ in P
    for ξ₂ in P
      if !(ϕ_2D(ξ₁, ξ₂) == ϕ_geral(ξ₁, ξ₂))
        teste_2D = false
        break
      end
    end
  end
  println(teste_2D ? "ok\n" : "nok\n")

end
# test_ϕ()
