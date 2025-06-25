function test_montaLG()
  ne = 5

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, ne)
  LG_1D = montaLG(ne, base)
  println("Teste montagem LG:")

  LG1D = LG_1D == montaLG_geral(5)
  println("LG 1D: $LG1D")

  LG2D = monta_LG_2D(5, 4) == montaLG_geral(5, 4)
  return println("LG 2D: $LG2D")
end
# test_montaLG()

function test_montaEQ()
  ne = 5

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, ne)
  neq = base.neq
  EQ_1D = montaEQ(ne, neq, base)
  println("\nTeste montagem EQ:")

  EQ1D = EQ_1D == montaEQ_geral(5)[2]
  println("EQ 1D: $EQ1D")
  # display(EQ_1D); display(montaEQ_geral(5))

  EQ2D = true
  # display(monta_EQ_2D(5, 4)[2] == montaEQ_geral(5, 4)[2])
  for nx1 in 1:50
    for nx2 in 1:50
      if !(monta_EQ_2D(nx1, nx2)[2] == montaEQ_geral(nx1, nx2)[2])
        EQ2D = false
      end
    end
  end
  return println("EQ 2D: $EQ2D")
end
# test_montaEQ()

function test_ϕ()
  println("\nTeste da ϕ:")
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
      # println("$(ϕ_2D(ξ₁, ξ₂)), $(ϕ_geral(ξ₁, ξ₂))")
      if !(ϕ_2D(ξ₁, ξ₂) == ϕ_geral(ξ₁, ξ₂)[1])
        teste_2D = false
        break
      end
    end
  end
  return println(teste_2D ? "ok\n" : "nok\n")
end
# test_ϕ()

function test_∇ϕ()
  println("Teste da ∇ϕ")
  P, W = legendre(2)

  print("∇ϕ 1D: ")
  teste_1D = true

  for ξ₁ in P
    if !(∇ϕ_1D(ξ₁) == ∇ϕ_geral(ξ₁))
      teste_1D = false
      break
    end
  end
  println(teste_1D ? "ok" : "nok")

  print("∇ϕ 2D: ")
  teste_2D = true

  for ξ₁ in P
    for ξ₂ in P
      ∂ϕ_∂ξ₁_2D, ∂ϕ_∂ξ₂_2D = ∇ϕ_2D(ξ₁, ξ₂)
      ∂ϕ_∂ξ₁_geral, ∂ϕ_∂ξ₂_geral = ∇ϕ_geral(ξ₁, ξ₂)
      if !(∂ϕ_∂ξ₁_2D == ∂ϕ_∂ξ₁_geral || ∂ϕ_∂ξ₂_2D == ∂ϕ_∂ξ₂_geral)
        teste_2D = false
        break
      end
    end
  end
  return println(teste_2D ? "ok\n" : "nok\n")
end
# test_∇ϕ()

function test_avalia_quadratura_∇ϕ()
  function avaliar_quadratura_geral_ϕ(
      base_func::Function, npg::Int64, n_funcs::Int64, n_dim::Int64)
    P, W = quadratura_gauss(npg, n_dim)

    ϕPₙ = ()
    for d in 1:n_dim
      ϕP = zeros(npg^n_dim, n_funcs^n_dim)
      for ξ in 1:(npg^n_dim)
        ϕP[ξ, :] .= base_func(P[ξ]...)[1]
      end
      ϕPₙ = (ϕPₙ..., ϕP)
    end

    return ϕPₙ, P, W
  end

  function avaliar_quadratura_geral(
      base_func::Function, npg::Int64, n_funcs::Int64, n_dim::Int64)
    P, W = quadratura_gauss(npg, n_dim)

    ϕPₙ = ()
    for d in 1:n_dim
      ϕP = zeros(npg^n_dim, n_funcs^n_dim)
      for ξ in 1:(npg^n_dim)
        ϕP[ξ, :] .= base_func(P[ξ]...)[d]
      end
      ϕPₙ = (ϕPₙ..., ϕP)
    end

    return ϕPₙ, P, W
  end

  println("Teste da avalia_quadratura para ∇ϕ")
  npg = 5

  print("ϕ 1D: ")
  n_dim = 1
  ∇ϕξ_1D, P_1D, W_1D = avaliar_quadratura_geral(∇ϕ_1D, npg, 2, n_dim)
  ∇ϕξ_geral, P_geral, W_geral = avaliar_quadratura_geral(∇ϕ_geral, npg, 2, n_dim)

  # display(∇ϕξ_1D)
  # display(∇ϕξ_geral)
  display(∇ϕξ_1D == ∇ϕξ_geral)

  print("ϕ 2D: ")
  n_dim = 2
  ∇ϕξ_2D, P_2D, W_2D = avaliar_quadratura_geral(∇ϕ_2D, npg, 2, n_dim)
  ∇ϕξ_geral, P_geral, W_geral = avaliar_quadratura_geral(∇ϕ_geral, npg, 2, n_dim)

  # display(∇ϕξ_2D)
  # display(∇ϕξ_geral)
  return display(∇ϕξ_2D == ∇ϕξ_geral)
  # display(P_geral)
end

# test_avalia_quadratura_∇ϕ()
