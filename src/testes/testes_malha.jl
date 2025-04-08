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
  display(EQ_1D == montaEQ_geral(5)[2])
  # display(EQ_1D); display(montaEQ_geral(5))

  println("\nEQ 2D")
  display(monta_EQ_2D(5, 4)[2] == montaEQ_geral(5, 4)[2])
  # display(monta_EQ_2D(5, 4)); display(montaEQ_geral(5, 4))
  for nx1 in 1:50
    for nx2 in 1:50
      if !(monta_EQ_2D(nx1, nx2)[2] == montaEQ_geral(nx1, nx2)[2])
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