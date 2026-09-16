using Fease
using LinearAlgebra
using GaussQuadrature
using SparseArrays

# ============== Início script =====================

# Infos Elementos Finitos
Nx1 = 20
Nx2 = 20

# Infos ALI
λ₁ = 1.0
λ₂ = 1.0

n_passos = 8
Δτ = 0.5 / n_passos

s1 = 1.0
s2 = -0.1
β = 10.0^4
params = [s1, s2, β]

# Infos Display
winX = 2.0
winY = 2.0

# Do domínio do problema
ponto_inf_esq = [0.0, 0.0]
ponto_sup_dir = [λ₁, λ₂]

# Neste exemplo, do quadrilatero apoiado no chão, estes são os nós prescritos
presc_x1 = [1]
presc_x2 = getQuadSideNodes(Nx1, Nx2, 1)
nos_prescritos = [presc_x1, presc_x2]

fronteira = monta_fronteira_2D_uniforme(
  ponto_inf_esq, ponto_sup_dir, Nx1, Nx2, nos_prescritos)

baseType = BaseTypes.linearLagrange
n_graus_liberdade = 2
malha = monta_malha_2D_uniforme(
  baseType, Nx1, Nx2, n_graus_liberdade, fronteira
)

P, W = legendre(2)
τ = 1 * Δτ
f = build_f(malha, P, τ)

# Construir quais valores estão prescritos. Ideal incluir estrutura de prescrição no objeto Malha
nodes = (unique(Iterators.flatten(fronteira.nos_prescritos)))
values = repeat([[0.0; 0.0]], size(nodes)[1])

g = build_g(nodes, malha.coords, values)

h = -s2 / s1
p = 1 - h
params = [s1, s2, p, β]

# F e T iniciais: Identidade e Matriz nula
F_def = repeat([Matrix(1.0I, 2, 2)], malha.ne)
T = repeat([Matrix(0.0I, 2, 2)], malha.ne)

X_novo = []
K_vec = []
F_vec = []
sol_vec = []
infos_primarias = []
infos_secundarias = []

for iter in 1:n_passos
  K, F = monta_K_F_global(params, f, g, malha, F_def, T)
  c = K \ F

  # u em t+1
  global X_novo = soma_solução(malha, c)
  #display(drawGrid(X_novo[:, 1], X_novo[:, 2], Nx1, Nx2, LG, winX, winY))

  novo_F, novo_T = atualiza_ALI(malha, X_novo, F_def, T, params)

  global F_def = copy(novo_F)
  global T = copy(novo_T)
  X = map(copy, X_novo)
  append!(sol_vec, [X])

  infos1,
  infos2 = calcula_erro(iter, F_def, s1, s2, τ, malha.fronteira.coords_quinas, malha.ne)
  append!(infos_primarias, [infos1])
  append!(infos_secundarias, [infos2])

  # Fim do passo atual, prepara pra novo passo

  novas_coords_quinas = [
    [ X[1][first(malha.fronteira.nos_fronteiras[1])],
      X[2][first(malha.fronteira.nos_fronteiras[1])]],

    [ X[1][last(malha.fronteira.nos_fronteiras[1])], 
      X[2][last(malha.fronteira.nos_fronteiras[1])]],

    [ X[1][last(malha.fronteira.nos_fronteiras[3])], 
      X[2][last(malha.fronteira.nos_fronteiras[3])]],

    [ X[1][first(malha.fronteira.nos_fronteiras[3])],
      X[2][first(malha.fronteira.nos_fronteiras[3])]]
  ]::Vector{Vector{Float64}}

  global fronteira = Fronteira(
    novas_coords_quinas,
    malha.fronteira.nos_fronteiras,
    malha.fronteira.elementos_fronteiras,
    malha.fronteira.nos_prescritos
  )

  global malha = Malha(
    malha.base,
    malha.ne, malha.neq,
    Tuple(X), malha.dx,
    malha.EQ, malha.LG, malha.EQoLG,
    malha.a, malha.b,
    malha.n_dim,
    malha.Nx,
    fronteira
  )

  global τ = Δτ * (iter + 1)
  global f = build_f(malha, P, τ)
end

#u = last(sol_vec)
#display(drawGrid(u[1], u[2], Nx1, Nx2, malha.LG, 1.8, 1.8))