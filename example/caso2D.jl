using Fease

Nx1, Nx2 = 4, 3
ne = Nx1 * Nx2

baseType = BaseTypes.linearLagrange
base = monta_base(baseType, ne)

a = (0.0, 0.0)
b = (1.0, 1.0)
malha = monta_malha_2D_uniforme(base, Nx1, Nx2, a, b)

example = 1
run_values = examples_2D(example)

(; α, β, f) = run_values

function pseudo_a(termos_equacao)
  (; ∇u, ∇v, u, v) = termos_equacao

  return β * dot(u, v) + α * dot(∇u, ∇v)
end

C = solve_sys(f, malha, pseudo_a)

display(C)

function convergence_2D()
  baseType = BaseTypes.linearLagrange

  a = (0.0, 0.0)
  b = (1.0, 1.0)

  example = 1
  run_values = examples_2D(example)

  errsize = 7
  NE = 2 .^ [2:1:errsize;]
  H = 1 ./ NE
  E = zeros(length(NE))
  dE = similar(E)

  convergence_test!(E, NE, run_values, baseType, a, b, 2)

  H = 1 ./ NE
  plot(H, E, xaxis = :log10, yaxis = :log10)
  return plot!(H, H .^ monta_base(baseType, 2).nB, xaxis = :log10, yaxis = :log10)
end

convergence_2D()
