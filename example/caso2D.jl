using MyProject

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
