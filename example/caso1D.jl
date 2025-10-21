using Fease

ne = 2^3

baseType = BaseTypes.linearLagrange
base = monta_base(baseType, ne)

a, b = 0, 1

malha = monta_malha_1D_uniforme(base, ne, a, b)

example = 1
run_values = examples_1D(example)

(; α, β, f) = run_values

function pseudo_a(termos_equacao)
  (; ∇u, ∇v, u, v) = termos_equacao

  return β * dot(u, v) + α * dot(∇u, ∇v)
end

C = solve_sys(f, malha, pseudo_a)

display(C)

function convergence_1D()
  baseType = BaseTypes.linearLagrange

  a, b = 0, 1

  example = 1
  run_values = examples_1D(example)

  errsize = 5
  NE = 2 .^ [2:1:errsize;]
  H = 1 ./ NE
  E = zeros(length(NE))
  dE = similar(E)
  convergence_test!(E, NE, run_values, baseType, a, b, 1)

  H = 1 ./ NE
  plot(H, E, xaxis = :log10, yaxis = :log10)
  return plot!(H, H .^ monta_base(baseType, 2).nB, xaxis = :log10, yaxis = :log10)
end

convergence_1D()