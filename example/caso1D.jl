using Fease

function caso_1D()
  # Define o número de sub-intervalos no eixo horizontal e vertical, respectivamente.
  Nx1 = 2^3

  # Define o tipo da base de funções interpoladoras do subespaço aproximado Vₘ.
  baseType = BaseTypes.linearLagrange

  # Define extremos do intervalo da malha.
  a, b = 0, 1

  # Define a malha com os valores atribuídos acima.
  malha = monta_malha_1D_uniforme(baseType, Nx1, a, b)

  # Define os parâmetros da equação a ser resolvida.
  example = 1
  run_values = examples_1D(example)
  (; α, β, f) = run_values

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
    return monta_malha_1D_uniforme(baseType, NX[1], a, b)
  end

  # Estudo será realizado com número de elementos em potências de 2.
  errsize = 7
  NE = 2 .^ [2:1:errsize;]
  H = 1 ./ NE
  E = zeros(length(NE))

  # Define a solução analítica da equação.
  (; u) = run_values

  # Calcula o erro para cada uma das quantidades de elementos finitos.
  n_dim = 1
  convergence_test!(E, NE, n_dim, monta_malha, pseudo_a, f, u)

  # Plota o resultado do estudo
  plot(H, E, xaxis = :log10, yaxis = :log10, label = "Erro")
  return plot!(H, H .^ monta_base(baseType, 2).nB, xaxis = :log10,
    yaxis = :log10, label = "H²")
end

caso_1D()