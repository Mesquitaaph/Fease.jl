# ========== Próprio do ALI ================

function calcL(s1::Float64, s2::Float64, β::Float64, B, inv_B)
  L = zeros(2, 2, 2, 2)
  delta = 1.0I

  for i in 1:2
    for j in 1:2
      for k in 1:2
        for l in 1:2
          L[i,
            j,
            k,
            l] = β * delta[k, l] * delta[i, j] +
                 s1 * (delta[i, k] * B[l, j] + B[i, l] * delta[j, k]) -
                 s2 * (inv_B[i, k] * delta[l, j] + delta[l, i] * inv_B[k, j])
        end
      end
    end
  end

  return L
end
export calcL

function calc_F_def_0(malha)
  F_def_0 = []

  ξ₁ = 0.0
  ξ₂ = 0.0
  ∇ϕξ = ∇ϕ_geral(ξ₁, ξ₂)

  for e in 1:malha.ne
    Xᵉ = elem_coords(malha, e)[2]

    dx_dξ₁ = mudanca_variavel_xξ(Xᵉ, ∇ϕξ[1], 2)
    dx_dξ₂ = mudanca_variavel_xξ(Xᵉ, ∇ϕξ[2], 2)

    F_defᵉ = [dx_dξ₁[1] dx_dξ₂[1]
              dx_dξ₁[2] dx_dξ₂[2]]

    append!(F_def_0, [F_defᵉ])
  end

  return F_def_0
end
export calc_F_def_0

function grad_u(malha_deformação, F_def_0)
  H = []

  for e in 1:malha_deformação.ne
    Δuᵉ = elem_coords(malha_deformação, e)[2]

    ξ₁ = 0.0
    ξ₂ = 0.0

    # quadratura_∇ϕ mas com P = vetor nulo e W = vetor unitário
    # Como é sempre o mesmo ponto, não precisaria calcular 4 vezes.
    # Mas estou fazendo esta estrutura pra reaproveitar o método de avaliação da ∇ϕ_geral
    P, W = [repeat([(0.0, 0.0)], 4), repeat([(1.0, 1.0)], 4)]
    n_funcs = malha_deformação.base.nB
    n_dim = malha_deformação.n_dim
    npg = 2

    ∇ϕP = ()
    for d in 1:n_dim
      ∂ϕᵢP = zeros(npg^n_dim, n_funcs^n_dim)

      # Para todos os pontos de Gauss, avalia as ∂ϕᵢ locais (i = d)
      for ξ in 1:(npg^n_dim)
        ∂ϕᵢP[ξ, :] .= ∇ϕ_geral(P[ξ]...)[d]
      end
      ∇ϕP = (∇ϕP..., ∂ϕᵢP)
    end

    # Tanto faz a linha, já que são todas iguais
    du_dξ₁ = mudanca_variavel_xξ(Δuᵉ, ∇ϕP[1][1, :], 2)
    du_dξ₂ = mudanca_variavel_xξ(Δuᵉ, ∇ϕP[2][1, :], 2)

    # reutilizando a função ∂ξ_to_∂x, mas idealmente o nome seria ∂ξ_to_∂u... ou "interpolate <something>"
    du_dξ = [du_dξ₁[1] du_dξ₂[1]
             du_dξ₁[2] du_dξ₂[2]]

    #display(du_dξ)
    dx_dξ = F_def_0[e]
    det_J_ξ = abs(det(dx_dξ))
    dξ_dx = inv(dx_dξ)

    Hᵉ = du_dξ * dξ_dx

    append!(H, [Hᵉ])
    #display(H[e])
  end

  return H
end
export grad_u

function atualiza_ALI(malha, X_novo, F_def, T, params)
  s1, s2, p, β = params     
  F_t_mais_1 = []
  T_t_mais_1 = []

  Δu = X_novo .- malha.coords
  malha_deformação = Malha(
    malha.base,
    malha.ne, malha.neq,
    Tuple(Δu), malha.dx,
    malha.EQ, malha.LG, malha.EQoLG,
    malha.a, malha.b,
    malha.n_dim,
    malha.Nx,
    malha.fronteira
  )

  F_def_0 = calc_F_def_0(malha)
  H = grad_u(malha_deformação, F_def_0)

  for e in 1:malha.ne
    F = F_def[e]

    det_J_ξ = abs(det(F))
    dξ_dx = inv(F)

    Hᵉ = H[e]
    B = F * transpose(F)
    inv_B = inv(B)

    L = calcL(s1, s2, β, B, inv_B)
    Tᵉ = copy(T[e])

    for i in 1:2
      for j in 1:2
        soma = 0.0
        for k in 1:2
          for l in 1:2
            Tᵉ[i, j] += L[i, j, k, l] * Hᵉ[k, l]
          end
        end
      end
    end

    # pressão
    #P_0 -= β * tr(H)

    Fᵉ_t_mais_1 = (1.0I + Hᵉ) * F

    #display(Fᵉ_t_mais_1)
    #display(Tᵉ)

    append!(F_t_mais_1, [Fᵉ_t_mais_1])
    append!(T_t_mais_1, [Tᵉ])
  end

  return F_t_mais_1, T_t_mais_1
end
export atualiza_ALI

# ============= Podem (e devem) ser mais gerais ==================

function bound_expr(num_fronteira, malha, P)
  nos_fronteira = malha.fronteira.nos_fronteiras[num_fronteira]
  e_fronteira = malha.fronteira.elementos_fronteiras[num_fronteira]
  X = malha.coords

  coords_nos_front = []
  for no in nos_fronteira
    append!(coords_nos_front, [[X[1][no], X[2][no]]])
  end

  for e in e_fronteira
    Xᵉ = elem_coords(malha, e)[2]

    if (num_fronteira == 1)
      for ξ in P
        x1_ξ1 = Xᵉ[1] ⋅ ϕ_2D(ξ, -1.0)
        x2_ξ2 = Xᵉ[2] ⋅ ϕ_2D(ξ, -1.0)
        append!(coords_nos_front, [[x1_ξ1, x2_ξ2]])
      end
    end

    if (num_fronteira == 2)
      for ξ in P
        x1_ξ1 = Xᵉ[1] ⋅ ϕ_2D(1.0, ξ)
        x2_ξ2 = Xᵉ[2] ⋅ ϕ_2D(1.0, ξ)
        append!(coords_nos_front, [[x1_ξ1, x2_ξ2]])
      end
    end

    if (num_fronteira == 3)
      for ξ in P
        x1_ξ1 = Xᵉ[1] ⋅ ϕ_2D(ξ, 1.0)
        x2_ξ2 = Xᵉ[2] ⋅ ϕ_2D(ξ, 1.0)
        append!(coords_nos_front, [[x1_ξ1, x2_ξ2]])
      end
    end

    if (num_fronteira == 4)
      for ξ in P
        x1_ξ1 = Xᵉ[1] ⋅ ϕ_2D(-1.0, ξ)
        x2_ξ2 = Xᵉ[2] ⋅ ϕ_2D(-1.0, ξ)
        append!(coords_nos_front, [[x1_ξ1, x2_ξ2]])
      end
    end
  end

  #println(coords_nos_front)

  expr = function (coords_no)
    if (coords_no in coords_nos_front)
      return true
    else
      return false
    end
  end
end
export bound_expr

function build_f(malha, P, τ::Float64)
  coords_quinas = malha.fronteira.coords_quinas
  nos_fronteira = malha.fronteira.nos_fronteiras
  e_fronteira = malha.fronteira.elementos_fronteiras

  p1 = coords_quinas[1]
  p2 = coords_quinas[2]
  p3 = coords_quinas[3]
  p4 = coords_quinas[4]

  dy = p4[2] - p1[2]
  dx = p4[1] - p1[1]
  dl = sqrt(dx^2 + dy^2)
  _sen = dy / dl
  _cos = dx / dl

  bound1 = bound_expr(1, malha, P)
  bound2 = bound_expr(2, malha, P)
  bound3 = bound_expr(3, malha, P)
  bound4 = bound_expr(4, malha, P)

  f = function (ponto)
    tracao = zeros(2)
    #println("Checking for ", ponto)

    if (bound1(ponto))
      tracao[1] += -τ
      tracao[2] += 0.0
    end

    if (bound2(ponto))
      tracao[1] += τ * _cos
      tracao[2] += τ * _sen
    end

    if (bound3(ponto))
      tracao[1] += τ
      tracao[2] += 0.0
    end

    if (bound4(ponto))
      tracao[1] += -τ * _cos
      tracao[2] += -τ * _sen
    end

    return tracao
  end

  return f
end
export build_f

function build_g(nodes, coords, values)::Function
  nodes_coords = []

  for node in nodes
    append!(nodes_coords, [coords[1][node], coords[2][node]])
  end

  g = function (x)
    if (x in nodes_coords)
      return values[findfirst(item -> item == x, nodes_coords)]
    else
      return [0.0; 0.0]
    end
  end

  return g
end
export build_g

function dirichlet_map(presc, e, LG)
  ndir = size(presc)[1]
  map = (zeros(Bool, 4), zeros(Bool, 4))
  for dir in 1:ndir
    for pt in 1:4
      if (LG[pt, e] in presc[dir])
        map[dir][pt] = true
      end
    end
  end

  return map
end
export dirichlet_map

# =========== Essas provavelmente serão extintas, para usar as do pacote =========

function quadratura_K_local(a::Int64, b::Int64, params::Array{Float64}, indexes, Xᵉ,
    F_defᵉ, Tᵉ, ∇ϕξ, P, W)::Float64
  quadratura = 0.0
  s1, s2, p, β = params
  r, s = indexes

  soma1 = 0.0
  soma2 = 0.0
  soma3 = 0.0

  n_combs = size(P)[1]
  for iter in 1:n_combs
    dx_dξ₁ = mudanca_variavel_xξ(Xᵉ, ∇ϕξ[1][iter, :], 2)
    dx_dξ₂ = mudanca_variavel_xξ(Xᵉ, ∇ϕξ[2][iter, :], 2)

    dx_dξ = [dx_dξ₁[1] dx_dξ₂[1]
             dx_dξ₁[2] dx_dξ₂[2]]

    det_J_ξ = abs(det(dx_dξ))
    dξ_dx = inv(dx_dξ)

    F = F_defᵉ
    B = F * transpose(F)
    inv_B = inv(B)

    L = calcL(s1, s2, β, B, inv_B)

    termo1 = 0.0
    termo2 = 0.0
    termo3 = 0.0

    ∇ϕ = ∇ϕ_2D(P[iter]...)

    for i in 1:2
      for k in 1:2
        for l in 1:2
          termo1 += Tᵉ[r, i] * ∇ϕ[k][a] * dξ_dx[k, i] * ∇ϕ[l][b] * dξ_dx[l, s] * det_J_ξ
          termo2 += Tᵉ[r, i] * ∇ϕ[k][a] * dξ_dx[k, s] * ∇ϕ[l][b] * dξ_dx[l, i] * det_J_ξ

          for j in 1:2
            termo3 += L[r, i, s, j] * ∇ϕ[k][a] * dξ_dx[k, i] * ∇ϕ[l][b] * dξ_dx[l, j] *
                      det_J_ξ
          end
        end
      end
    end

    soma1 += termo1
    soma2 += termo2
    soma3 += termo3

    quadratura += reduce(*, W[iter]) * (termo1 - termo2 + termo3)
  end

  return quadratura
end

function monta_K_local(params::Array{Float64}, Xᵉ, F_defᵉ, Tᵉ,
    ∇ϕξ, P, W)::Matrix{Float64}
  Kᵉ = zeros(8, 8)

  for a in 1:4
    for r in 1:2
      for b in 1:4
        for s in 1:2
          m, n = (2 * (a - 1) + r, 2 * (b - 1) + s)

          indexes = [r, s]
          Kᵉ[m, n] += quadratura_K_local(a, b, params, indexes, Xᵉ, F_defᵉ, Tᵉ, ∇ϕξ, P, W)
        end
      end
    end
  end

  return Kᵉ
end

function quadratura_F_local(
    f::Function, g::Function, presc_map, a::Int64, indexes, boundaries, Xᵉ,
    Kᵉ::Matrix{Float64}, F_defᵉ, Tᵉ, ∇ϕξ, P::Vector{Float64}, W::Vector{Float64}, paired_P, paired_W
)
  quadratura = 0.0
  r = indexes

  if (length(boundaries) != 0)
    for b in 1:4
      for (ξ, w) in zip(P, W)
        termo1 = 0.0

        # Integral lado 1
        p1 = [Xᵉ[1][1], Xᵉ[2][1]]
        p2 = [Xᵉ[1][2], Xᵉ[2][2]]

        x1_ξ1 = Xᵉ[1] ⋅ ϕ_2D(ξ, -1.0)
        x2_ξ2 = Xᵉ[2] ⋅ ϕ_2D(ξ, -1.0)

        tamanho_lado = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        ϕ = ϕ_2D(ξ, -1.0)

        termo1 += w * ϕ[a] * ϕ[b] * f([x1_ξ1, x2_ξ2])[r] * tamanho_lado/2
        # fator 2 = tamanho do lado 
        # do elemento padrão

        # Integral lado 2
        p1 = [Xᵉ[1][2], Xᵉ[2][2]]
        p2 = [Xᵉ[1][3], Xᵉ[2][3]]

        x1_ξ1 = Xᵉ[1] ⋅ ϕ_2D(1.0, ξ)
        x2_ξ2 = Xᵉ[2] ⋅ ϕ_2D(1.0, ξ)

        tamanho_lado = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        ϕ = ϕ_2D(1.0, ξ)

        termo1 += w * ϕ[a] * ϕ[b] * f([x1_ξ1, x2_ξ2])[r] * tamanho_lado/2

        # Integral lado 3
        p1 = [Xᵉ[1][3], Xᵉ[2][3]]
        p2 = [Xᵉ[1][4], Xᵉ[2][4]]

        x1_ξ1 = Xᵉ[1] ⋅ ϕ_2D(ξ, 1.0)
        x2_ξ2 = Xᵉ[2] ⋅ ϕ_2D(ξ, 1.0)

        tamanho_lado = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        ϕ = ϕ_2D(ξ, 1.0)

        termo1 += w * ϕ[a] * ϕ[b] * f([x1_ξ1, x2_ξ2])[r] * tamanho_lado/2

        # Integral lado 4
        p1 = [Xᵉ[1][1], Xᵉ[2][1]]
        p2 = [Xᵉ[1][4], Xᵉ[2][4]]

        x1_ξ1 = Xᵉ[1] ⋅ ϕ_2D(-1.0, ξ)
        x2_ξ2 = Xᵉ[2] ⋅ ϕ_2D(-1.0, ξ)

        tamanho_lado = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        ϕ = ϕ_2D(-1.0, ξ)

        termo1 += w * ϕ[a] * ϕ[b] * f([x1_ξ1, x2_ξ2])[r] * tamanho_lado/2

        quadratura += termo1
      end
    end
  end

  prescrito = presc_map[r][a]

  if (!prescrito)
    n_combs = size(paired_P)[1]
    for i in 1:2
      soma = 0.0
      somaP = 0.0
      for iter in 1:n_combs
        dx_dξ₁ = mudanca_variavel_xξ(Xᵉ, ∇ϕξ[1][iter, :], 2)
        dx_dξ₂ = mudanca_variavel_xξ(Xᵉ, ∇ϕξ[2][iter, :], 2)

        dx_dξ = [dx_dξ₁[1] dx_dξ₂[1]
                 dx_dξ₁[2] dx_dξ₂[2]]

        det_J_ξ = abs(det(dx_dξ))
        dξ_dx = inv(dx_dξ)

        ∇ϕ = ∇ϕ_2D(paired_P[iter]...)
        termo2 = 0.0
        for k in 1:2
          termo2 += reduce(*, paired_W[iter]) * Tᵉ[r, i] * ∇ϕ[k][a] * dξ_dx[k, i] * det_J_ξ
          somaP += ∇ϕ[k][a] * dξ_dx[k, i] * det_J_ξ
        end

        soma += termo2
        quadratura -= termo2
      end
    end
  end

  for b in 1:4
    ponto = (Xᵉ[1][b], Xᵉ[2][b])

    termo3 = 0.0
    termo3 += Kᵉ[a, b] * g(ponto)[r]

    quadratura -= termo3
  end

  return quadratura
end

function monta_F_local(f::Function, g, presc_map, boundaries, Xᵉ, Kᵉ, F_defᵉ, Tᵉ, ∇ϕξ,
    P::Vector{Float64}, W::Vector{Float64}, paired_P, paired_W)::Vector{Float64}
  F_local = zeros(8)

  for a in 1:4
    for r in 1:2
      m = 2 * (a - 1) + r

      indexes = r
      F_local[m] += quadratura_F_local(
        f, g, presc_map, a, indexes, boundaries, Xᵉ, Kᵉ,
        F_defᵉ, Tᵉ, ∇ϕξ, P, W, paired_P, paired_W)
    end
  end

  return F_local
end


function monta_K_F_global(params::Array{Float64}, f::Function, g::Function, malha, F_def, T)
  m = malha.neq

  K_global = spzeros(m + 1, m + 1)
  F_global = zeros(m + 1)

  # neq dessa estrutura só vale quando não tem prescrição, usar neq da malha como na primeira linha
  base = monta_base(BaseTypes.linearLagrange, malha.ne)

  ∇ϕξ, paired_P, paired_W = quadratura_∇ϕ(base, 2, 2)
  P, W = legendre(2)

  for e in 1:malha.ne
    #println("e = ", e)
    e_boundaries = getBoundariesOfElement(e, malha.Nx[1], malha.Nx[2], malha.LG)

    Xᵉ = elem_coords(malha, e)[2]

    F_defᵉ = F_def[e]
    Tᵉ = T[e]
    # Tirar essas 2 linhas pra fora do for ao final, depois de ajustar a monta_F
    Kᵉ = monta_K_local(params, Xᵉ, F_defᵉ, Tᵉ, ∇ϕξ, paired_P, paired_W)

    presc_map = dirichlet_map(malha.fronteira.nos_prescritos, e, malha.LG)

    Fᵉ = monta_F_local(
      f, g, presc_map, e_boundaries, Xᵉ, Kᵉ, F_defᵉ, Tᵉ, ∇ϕξ, P, W, paired_P, paired_W)

    for a in 1:4
      for r in 1:2
        index_i = malha.EQ[:, r][malha.LG][a, e]

        for b in 1:4
          for s in 1:2
            index_j = malha.EQ[:, s][malha.LG][b, e]
            K_global[index_i, index_j] += Kᵉ[2*(a-1)+r, 2*(b-1)+s]
          end
        end
        F_global[index_i] += Fᵉ[2*(a-1)+r]
      end
    end
  end

  return (K_global[1:m, 1:m], F_global[1:m])
end
export monta_K_F_global

# Não é a mesma coisa que a monta_u_aproximada!
function soma_solução(malha, solução)
  X_novo = map(copy, malha.coords)
  npts = (malha.Nx[1] + 1) * (malha.Nx[2] + 1)

  for pt in 1:npts
    for dir in 1:2
      idx = malha.EQ[pt, dir]
      if idx != malha.neq + 1
        X_novo[dir][pt] += solução[idx]
      end
    end
  end

  return X_novo
end
export soma_solução