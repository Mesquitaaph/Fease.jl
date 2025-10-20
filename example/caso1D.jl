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
