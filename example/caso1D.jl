using Fease

ne = 2^3

baseType = BaseTypes.linearLagrange

a, b = 0, 1

malha = monta_malha_1D_uniforme(baseType, ne, a, b)

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

  function monta_malha(NX)
    return monta_malha_1D_uniforme(baseType, NX[1], a, b)
  end

  example = 1
  run_values = examples_1D(example)

  errsize = 5
  NE = 2 .^ [2:1:errsize;]
  H = 1 ./ NE
  E = zeros(length(NE))
  dE = similar(E)
  convergence_test!(E, NE, run_values, 1, monta_malha)

  H = 1 ./ NE
  plot(H, E, xaxis = :log10, yaxis = :log10)
  return plot!(H, H .^ monta_base(baseType, 2).nB, xaxis = :log10, yaxis = :log10)
end

convergence_1D()
