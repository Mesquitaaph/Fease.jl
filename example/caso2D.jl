using Fease

function caso_2D()
  # Define o número de sub-intervalos no eixo horizontal e vertical, respectivamente.
  Nx1, Nx2 = 4, 4

  # Define o tipo da base de funções interpoladoras do subespaço aproximado Vₘ.
  baseType = BaseTypes.linearLagrange

  # Define extremos do intervalo da malha.
  a = (0.0, 0.0)
  b = (1.0, 1.0)

  # Define a malha com os valores atribuídos acima.
  malha = monta_malha_2D_uniforme(baseType, Nx1, Nx2, a, b)

  # Define os parâmetros da equação a ser resolvida.
  α = 1.0
  β = 1.0
  f = (x₁, x₂) -> (2 * α * π^2 + β) * sin(π * x₁) * sin(π * x₂)

  # Define o pseudo operador linear a(u,v).
  function pseudo_a(termos_equacao)
    (; ∇u, ∇v, u, v) = termos_equacao

    return β * dot(u, v) + α * dot(∇u, ∇v)
  end

  # Monta e resolve o sistema linear relacionado a esta equação.
  c = solve_sys(f, malha, pseudo_a)

  # Imprime solução aproximada.
  show(c)

  # Plota solução aproximada.
  plot_solucao_aproximada(c, malha, false)

  # Estudo da convergencia do erro

  # Define função que constrói uma malha que depende só do número de sub-intervalos nos dois eixos.
  function monta_malha(NX)
    return monta_malha_2D_uniforme(baseType, NX..., a, b)
  end

  # Estudo será realizado com número de elementos em potências de 2.
  errsize = 7
  NE = 2 .^ [2:1:errsize;]
  H = 1 ./ NE
  E = zeros(length(NE))

  # Define a solução analítica da equação.
  u = (x₁, x₂) -> sin(π * x₁) * sin(π * x₂)

  # Calcula o erro para cada uma das quantidades de elementos finitos.
  convergence_test!(E, NE, 2, monta_malha, pseudo_a, f, u)

  # Plota o resultado do estudo
  plot(H, E, xaxis = :log10, yaxis = :log10, label = "Erro")
  return plot!(H, H .^ monta_base(baseType, 2).nB, xaxis = :log10,
    yaxis = :log10, label = "H²")
end

# Chama a função definida acima.
caso_2D()
