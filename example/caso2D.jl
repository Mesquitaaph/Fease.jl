using Fease

Nx1, Nx2 = 4, 4
ne = Nx1 * Nx2

baseType = BaseTypes.linearLagrange

a = (0.0, 0.0)
b = (1.0, 1.0)
malha = monta_malha_2D_uniforme(baseType, Nx1, Nx2, a, b)

alpha = 1.0
beta = 1.0
f = (x₁, x₂) -> (2 * alpha * π^2 + beta) * sin(π * x₁) * sin(π * x₂)
u = (x₁, x₂) -> sin(π * x₁) * sin(π * x₂)

function ref_op_a(termos_equacao)
  (; ∇u, ∇v, u, v) = termos_equacao

  return β * dot(u, v) + α * dot(∇u, ∇v)
end

C = solve_sys(f, malha, pseudo_a)

display(C)

plot_solucao_aproximada(C, malha, false)

function convergence_2D()
  baseType = BaseTypes.linearLagrange

  a = (0.0, 0.0)
  b = (1.0, 1.0)

  function monta_malha(NX)
    return monta_malha_2D_uniforme(baseType, NX..., a, b)
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

# convergence_2D()
