using Fease
using LinearAlgebra
using GaussQuadrature
using SparseArrays

# ============== Funções específicas ===============
function mapper_to_x_generic(Xᵉ_a::Vector{Float64})::Function
  return (ξ₁::Float64, ξ₂::Float64)->Xᵉ_a ⋅ map(f -> f(ξ₁, ξ₂), ϕ.(1:length(Xᵉ_a)))
end

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
    ξ_to_x1 = mapper_to_x_generic(Xᵉ[1])
    ξ_to_x2 = mapper_to_x_generic(Xᵉ[2])

    if (num_fronteira == 1)
      for ξ in P
        x1_ξ1 = ξ_to_x1(ξ, -1.0)
        x2_ξ2 = ξ_to_x2(ξ, -1.0)
        append!(coords_nos_front, [[x1_ξ1, x2_ξ2]])
      end
    end

    if (num_fronteira == 2)
      for ξ in P
        x1_ξ1 = ξ_to_x1(1.0, ξ)
        x2_ξ2 = ξ_to_x2(1.0, ξ)
        append!(coords_nos_front, [[x1_ξ1, x2_ξ2]])
      end
    end

    if (num_fronteira == 3)
      for ξ in P
        x1_ξ1 = ξ_to_x1(ξ, 1.0)
        x2_ξ2 = ξ_to_x2(ξ, 1.0)
        append!(coords_nos_front, [[x1_ξ1, x2_ξ2]])
      end
    end

    if (num_fronteira == 4)
      for ξ in P
        x1_ξ1 = ξ_to_x1(-1.0, ξ)
        x2_ξ2 = ξ_to_x2(-1.0, ξ)
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

function simpleDet(M)
  return (M[1, 1] * M[2, 2]) - (M[1, 2] * M[2, 1])
end

function simpleInv(M)
  det_M = abs(simpleDet(M))
  inv_M = 1 / det_M * [M[2, 2] -M[1, 2]
                       -M[2, 1] M[1, 1]]
  return inv_M
end

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

function ϕ(a::Int64)::Function
  if a == 1
    return (ξ₁::Float64, ξ₂::Float64) -> (1.0-ξ₁)*(1.0-ξ₂)/4.0

  elseif a == 2
    return (ξ₁::Float64, ξ₂::Float64) -> (1.0+ξ₁)*(1.0-ξ₂)/4.0

  elseif a == 3
    return (ξ₁::Float64, ξ₂::Float64) -> (1.0+ξ₁)*(1.0+ξ₂)/4.0

  elseif a == 4
    return (ξ₁::Float64, ξ₂::Float64) -> (1.0-ξ₁)*(1.0+ξ₂)/4.0
  end
end

function ∂ϕ(a::Int64, variable::Int64)::Function
  if variable == 1
    if a == 1
      return (ξ₁::Float64, ξ₂::Float64) -> -(1.0-ξ₂)/4.0
    elseif a == 2
      return (ξ₁::Float64, ξ₂::Float64) -> (1.0-ξ₂)/4.0
    elseif a == 3
      return (ξ₁::Float64, ξ₂::Float64) -> (1.0+ξ₂)/4.0
    elseif a == 4
      return (ξ₁::Float64, ξ₂::Float64) -> -(1.0+ξ₂)/4.0
    end

  elseif variable == 2
    if a == 1
      return (ξ₁::Float64, ξ₂::Float64) -> -(1.0-ξ₁)/4.0
    elseif a == 2
      return (ξ₁::Float64, ξ₂::Float64) -> -(1.0+ξ₁)/4.0
    elseif a == 3
      return (ξ₁::Float64, ξ₂::Float64) -> (1.0+ξ₁)/4.0
    elseif a == 4
      return (ξ₁::Float64, ξ₂::Float64) -> (1.0-ξ₁)/4.0
    end
  end
end

function ξ_to_x(Xᵉ_a::Vector{Float64}, ξ₁::Float64, ξ₂::Float64)::Float64
  ϕ_ξ = map(func -> func(ξ₁, ξ₂), ϕ.(1:length(Xᵉ_a)))

  return (Xᵉ_a ⋅ ϕ_ξ)
end

function ∂ξ_to_∂x(Xᵉ_a::Vector{Float64}, variable::Int64, ξ₁::Float64, ξ₂::Float64)::Float64
  dϕ_dξ = map(func -> func(ξ₁, ξ₂), ∂ϕ.(1:length(Xᵉ_a), variable))

  return (Xᵉ_a ⋅ dϕ_dξ)
end

function elementCoords(element, LG, X)
  return [X[1][LG[:, element]];; X[2][LG[:, element]]]
end

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

    #display(B)
    #display(inv_B)

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
          #println("a, b = ", a, b, "\n Indexes (r, s, i, j) = ", indexes)
          Kᵉ[m, n] += quadratura_K_local(a, b, params, indexes, Xᵉ, F_defᵉ, Tᵉ, ∇ϕξ, P, W)
        end
      end
    end
  end

  return Kᵉ
end

function quadratura_F_local(
    f::Function, g::Function, presc_map, a::Int64, indexes, boundaries, Xᵉ,
    Kᵉ::Matrix{Float64}, F_defᵉ, Tᵉ, P::Vector{Float64}, W::Vector{Float64}
)
  quadratura = 0.0
  r = indexes

  ξ_to_x1 = mapper_to_x_generic(Xᵉ[1])
  ξ_to_x2 = mapper_to_x_generic(Xᵉ[2])

  if (length(boundaries) != 0)
    for b in 1:4
      ponto = (Xᵉ[1][b], Xᵉ[2][b])
      #display(f(ponto))

      for (ξ, w) in zip(P, W)
        termo1 = 0.0

        # Integral lado 1
        p1 = [Xᵉ[1][1], Xᵉ[2][1]]
        p2 = [Xᵉ[1][2], Xᵉ[2][2]]

        x1_ξ1 = ξ_to_x1(ξ, -1.0)
        x2_ξ2 = ξ_to_x2(ξ, -1.0)

        tamanho_lado = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        termo1 += w * ϕ(a)(ξ, -1.0) * ϕ(b)(ξ, -1.0) * f([x1_ξ1, x2_ξ2])[r] * tamanho_lado/2
        # fator 2 = tamanho do lado 
        # do elemento padrão
        #println("f1: ", f([x1_ξ1, x2_ξ2])[r])
        # Integral lado 2
        p1 = [Xᵉ[1][2], Xᵉ[2][2]]
        p2 = [Xᵉ[1][3], Xᵉ[2][3]]

        x1_ξ1 = ξ_to_x1(1.0, ξ)
        x2_ξ2 = ξ_to_x2(1.0, ξ)

        tamanho_lado = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        termo1 += w * ϕ(a)(1.0, ξ) * ϕ(b)(1.0, ξ) * f([x1_ξ1, x2_ξ2])[r] * tamanho_lado/2

        #println("f2: ", f([x1_ξ1, x2_ξ2])[r])
        # Integral lado 3
        p1 = [Xᵉ[1][3], Xᵉ[2][3]]
        p2 = [Xᵉ[1][4], Xᵉ[2][4]]

        x1_ξ1 = ξ_to_x1(ξ, 1.0)
        x2_ξ2 = ξ_to_x2(ξ, 1.0)

        tamanho_lado = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        termo1 += w * ϕ(a)(ξ, 1.0) * ϕ(b)(ξ, 1.0) * f([x1_ξ1, x2_ξ2])[r] * tamanho_lado/2

        #println("f3: ", f([x1_ξ1, x2_ξ2])[r])
        # Integral lado 4
        p1 = [Xᵉ[1][1], Xᵉ[2][1]]
        p2 = [Xᵉ[1][4], Xᵉ[2][4]]

        x1_ξ1 = ξ_to_x1(-1.0, ξ)
        x2_ξ2 = ξ_to_x2(-1.0, ξ)

        tamanho_lado = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        termo1 += w * ϕ(a)(-1.0, ξ) * ϕ(b)(-1.0, ξ) * f([x1_ξ1, x2_ξ2])[r] * tamanho_lado/2

        #println("f4: ", f([x1_ξ1, x2_ξ2])[r])
        quadratura += termo1
        #println("quad = ", quadratura, " a = ", a, " b = ", b)
      end
    end
  end

  prescrito = presc_map[r][a]

  if (!prescrito)
    n_pts = size(P)[1]
    for i in 1:2
      soma = 0.0
      somaP = 0.0
      for iter1 in 1:n_pts
        for iter2 in 1:n_pts
          ξ₁, w₁ = P[iter1], W[iter1]
          ξ₂, w₂ = P[iter2], W[iter2]

          dx_dξ = [∂ξ_to_∂x(Xᵉ[1], 1, ξ₁, ξ₂) ∂ξ_to_∂x(Xᵉ[1], 2, ξ₁, ξ₂)
                   ∂ξ_to_∂x(Xᵉ[2], 1, ξ₁, ξ₂) ∂ξ_to_∂x(Xᵉ[2], 2, ξ₁, ξ₂)]

          #display(dx_dξ)
          det_J_ξ = abs(simpleDet(dx_dξ))
          dξ_dx = simpleInv(dx_dξ)

          termo2 = 0.0
          for k in 1:2
            termo2 += (w₁ * w₂) * Tᵉ[r, i] * ∂ϕ(a, k)(ξ₁, ξ₂) * dξ_dx[k, i] * det_J_ξ
            # println("P; k = ", k, ": ", ∂ϕ(a, k)(ξ₁, ξ₂) * dξ_dx[k, i] * det_J_ξ)
            # println("dPhi = ", ∂ϕ(a, k)(ξ₁, ξ₂))
            # println("dξ_dx = ", dξ_dx[k, i])
            # println("det_J_ξ = ", det_J_ξ)
            # println("")
            somaP += ∂ϕ(a, k)(ξ₁, ξ₂) * dξ_dx[k, i] * det_J_ξ
          end
          soma += termo2
          quadratura -= termo2
        end
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

function monta_F_local(f::Function, g, presc_map, boundaries, Xᵉ, Kᵉ, F_defᵉ, Tᵉ,
    P::Vector{Float64}, W::Vector{Float64})::Vector{Float64}
  F_local = zeros(8)

  for a in 1:4
    for r in 1:2
      m = 2 * (a - 1) + r

      indexes = r
      F_local[m] += quadratura_F_local(
        f, g, presc_map, a, indexes, boundaries, Xᵉ, Kᵉ, F_defᵉ, Tᵉ, P, W)
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

  for e in 1:malha.ne
    #println("e = ", e)
    e_boundaries = getBoundariesOfElement(e, malha.Nx[1], malha.Nx[2], malha.LG)

    Xᵉ = elem_coords(malha, e)[2]

    F_defᵉ = F_def[e]
    Tᵉ = T[e]

    ϕξ, P, W = quadratura_ϕ(base, 2, 2)
    ∇ϕξ, P, W = quadratura_∇ϕ(base, 2, 2)
    # Tirar essas 2 linhas pra fora do for ao final, depois de ajustar a monta_F
    Kᵉ = monta_K_local(params, Xᵉ, F_defᵉ, Tᵉ, ∇ϕξ, P, W)

    presc_map = dirichlet_map(malha.fronteira.nos_prescritos, e, malha.LG)

    P, W = legendre(2) # enquanto não ajeitar a monta_F, usar padrão antigo
    Fᵉ = monta_F_local(f, g, presc_map, e_boundaries, Xᵉ, Kᵉ, F_defᵉ, Tᵉ, P, W)
    #display(Fᵉ)
    #println("\n", "Elemento ", e, "\n")
    for a in 1:4
      for r in 1:2
        index_i = malha.EQ[:, r][malha.LG][a, e]
        for b in 1:4
          for s in 1:2
            index_j = malha.EQ[:, s][malha.LG][b, e]
            #println("[", index_i, ", ", index_j, "] com ", "a = ", a, ", r = ", r, " b = ", b, ", s = ", s)
            K_global[index_i, index_j] += Kᵉ[2*(a-1)+r, 2*(b-1)+s]
          end
        end
        F_global[index_i] += Fᵉ[2*(a-1)+r]
        #println("F[", index_i, "] = ", F_global[index_i])
      end
    end
  end

  return (K_global[1:m, 1:m], F_global[1:m])
end

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

function calc_F_def_0(malha)
  F_def_0 = []

  ξ₁ = 0.0
  ξ₂ = 0.0
  for e in 1:malha.ne
    Xᵉ = elem_coords(malha, e)[2]
    F_defᵉ = [∂ξ_to_∂x(Xᵉ[1], 1, ξ₁, ξ₂) ∂ξ_to_∂x(Xᵉ[1], 2, ξ₁, ξ₂)
              ∂ξ_to_∂x(Xᵉ[2], 1, ξ₁, ξ₂) ∂ξ_to_∂x(Xᵉ[2], 2, ξ₁, ξ₂)]

    append!(F_def_0, [F_defᵉ])
  end

  return F_def_0
end

function grad_u(Δu, malha, F_def_0)
  H = []

  for e in 1:malha.ne
    Δuᵉ = elementCoords(e, malha.LG, Δu)

    ξ₁ = 0.0
    ξ₂ = 0.0

    # reutilizando a função ∂ξ_to_∂x, mas idealmente o nome seria ∂ξ_to_∂u... ou "interpolate <something>"
    du_dξ = [∂ξ_to_∂x(Δuᵉ[:, 1], 1, ξ₁, ξ₂) ∂ξ_to_∂x(Δuᵉ[:, 1], 2, ξ₁, ξ₂)
             ∂ξ_to_∂x(Δuᵉ[:, 2], 1, ξ₁, ξ₂) ∂ξ_to_∂x(Δuᵉ[:, 2], 2, ξ₁, ξ₂)]

    dx_dξ = F_def_0[e]

    #display("du_dξ = ")
    #display(du_dξ)
    det_J_ξ = abs(simpleDet(dx_dξ))
    dξ_dx = simpleInv(dx_dξ)

    Hᵉ = du_dξ * dξ_dx

    append!(H, [Hᵉ])
    #display(H[e])
  end

  return H
end

function atualiza_ALI(malha, X_novo, F_def, T)
  F_t_mais_1 = []
  T_t_mais_1 = []

  Δu = X_novo .- malha.coords

  F_def_0 = calc_F_def_0(malha)
  H = grad_u(Δu, malha, F_def_0)

  for e in 1:malha.ne
    F = F_def[e]

    det_J_ξ = abs(simpleDet(F))
    dξ_dx = simpleInv(F)

    Hᵉ = H[e]
    B = F * transpose(F)
    inv_B = simpleInv(B)

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

function calcula_erro(passo, F_def, s1, s2, τ, pts_fronteira, ne)
  max_F11 = 0.0
  min_F11 = 2.0
  sum_F11 = 0.0
  sum_det = 0.0

  for F in F_def
    sum_det += simpleDet(F)
    sum_F11 += F[1, 1]
    max_F11 = (F[1, 1] > max_F11) ? F[1, 1] : max_F11
    min_F11 = (F[1, 1] < min_F11) ? F[1, 1] : min_F11
  end

  λ₁_analitico = (1.0 - τ^2 * (1.0 / (s1 - s2)^2))^(-1.0 / 4)
  k = sqrt(λ₁_analitico^4 - 1.0)

  p1 = pts_fronteira[1]
  p2 = pts_fronteira[2]
  p4 = pts_fronteira[4]

  OA_num = p2[1] - p1[1]
  OB_num = sqrt((p4[1]^2 - p1[1]^2) + (p4[2]^2 - p1[2]^2))

  OA_analitico = λ₁_analitico^2
  #F_11 médio
  λ₁_numerico = sum_F11 / ne
  det_medio = sum_det / ne

  erro_relativo_percentual = 100 * abs((λ₁_numerico - λ₁_analitico) / λ₁_analitico)
  razao_segmentos = 100 * abs((OA_num - OB_num) / OA_num)
  infos_principais = [passo,
    λ₁_analitico,
    λ₁_numerico,
    min_F11,
    max_F11,
    det_medio,
    erro_relativo_percentual,
    razao_segmentos]::Vector{Float64}

  infos_secundarias = [k, OA_analitico]

  return infos_principais, infos_secundarias
end

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

  novo_F, novo_T = atualiza_ALI(malha, X_novo, F_def, T)

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
    [X[1][first(malha.fronteira.nos_fronteiras[1])],
      X[2][first(malha.fronteira.nos_fronteiras[1])]],
    [
      X[1][last(malha.fronteira.nos_fronteiras[1])], X[2][last(malha.fronteira.nos_fronteiras[1])]],
    [
      X[1][last(malha.fronteira.nos_fronteiras[3])], X[2][last(malha.fronteira.nos_fronteiras[3])]],
    [X[1][first(malha.fronteira.nos_fronteiras[3])],
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