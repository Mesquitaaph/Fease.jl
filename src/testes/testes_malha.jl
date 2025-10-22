function teste_foco_2D()
  baseType = BaseTypes.linearLagrange

  a = (0.0, 0.0)
  b = (1.0, 1.0)

  Nx1, Nx2 = 4, 3
  ponto_foco = (0.3, 0.2)

  # malha = monta_malha_2D_uniforme(baseType, Nx1, Nx2, a, b)
  malha = monta_malha_2D_foco(baseType, Nx1, Nx2, a, b, ponto_foco, 5)

  plot_malha_2D(malha)
  example = 1
  run_values = examples_2D(example)

  c = solve_sys_poisson(run_values, malha)

  # plot_solucao_aproximada(c, malha, false)

  return
end

teste_foco_2D()

function convergence_2D()
  baseType = BaseTypes.linearLagrange

  a = (0.0, 0.0)
  b = (1.0, 1.0)

  function monta_malha(NX)
    ponto_foco = (0.7, 0.5)
    return monta_malha_2D_foco(baseType, NX..., a, b, ponto_foco, 3)
    # return monta_malha_2D_uniforme(baseType, NX..., a, b)
  end

  example = 1
  run_values = examples_2D(example)

  errsize = 7
  NE = 2 .^ [2:1:errsize;]
  H = 1 ./ NE
  E = zeros(length(NE))
  dE = similar(E)

  convergence_test!(E, NE, run_values, 2, monta_malha)

  H = 1 ./ NE
  plot(H, E, xaxis = :log10, yaxis = :log10)
  return plot!(H, H .^ monta_base(baseType, 2).nB, xaxis = :log10, yaxis = :log10)
end

convergence_2D()
