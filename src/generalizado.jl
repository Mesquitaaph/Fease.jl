function jacobiano(n_dim, Xe, ∇Φξ, ξ)
  M = Matrix{Float64}(undef, n_dim, n_dim)
  for i in 1:n_dim
    for j in 1:n_dim
      M[i,j] = dot(Xe[i], ∇Φξ[j][ξ,:])
    end
  end
  detJ = det(M)

  return M, detJ
end

function quadratura_gauss(npg::Int64, n_dim::Int64)
  p, w = legendre(npg)

  P = Vector{Tuple}(undef, npg^n_dim)
  W = similar(P)

  for i in 1:npg^n_dim
    ponto = ()
    peso = ()
    for d in 1:n_dim
      idx_global = floor((Float64(npg)^(n_dim-d) + i-1) * Float64(npg)^-(n_dim-d))
      idx = Int((idx_global-1) % npg)+1
      ponto = (ponto..., p[idx])
      peso = (peso..., w[idx])
    end
    P[i] = ponto
    W[i] = peso
  end

  return P, W
end

function quadratura_ϕ(base, npg, n_dim)
  P, W = quadratura_gauss(npg, n_dim)
  n_funcs = base.nB

  ϕPₙ = ()
  for d in 1:n_dim
    ϕP = zeros(npg^n_dim, n_funcs^n_dim)
    for ξ in 1:npg^n_dim
      ϕP[ξ, :] .= ϕ_geral(P[ξ]...)[1]
    end
    ϕPₙ = (ϕPₙ..., ϕP)
  end

  return ϕPₙ, P, W
end

function quadratura_∇ϕ(base, npg, n_dim)
  P, W = quadratura_gauss(npg, n_dim)
  n_funcs = base.nB

  ϕPₙ = ()
  for d in 1:n_dim
    ϕP = zeros(npg^n_dim, n_funcs^n_dim)
    for ξ in 1:npg^n_dim
      ϕP[ξ, :] .= ∇ϕ_geral(P[ξ]...)[d]
    end
    ϕPₙ = (ϕPₙ..., ϕP)
  end

  return ϕPₙ, P, W
end

function montaKᵉ_geral!(Kᵉ, Xe, P, W, Φξ, ∇Φξ, n_dim, dx, run_values::RunValues)
  fill!(Kᵉ, 0.0)

  function ∇Φ(ξ, a)
    return [∇Φξ[dim][ξ, a] for dim in 1:n_dim]
  end
  
  (; α, β, γ) = run_values

  npg = length(P)
  for ξ in 1:npg
    ϕᵉ = Φξ[1][ξ,:]

    M, detJ = jacobiano(n_dim, Xe, ∇Φξ, ξ)
    @assert detJ > 0 "O determinante jacobiano deve ser positivo"
    
    M⁻¹ = inv(M)
    Hᵀ = M⁻¹*detJ
    H = transpose(Hᵀ)

    WW = prod(W[ξ])
    
    for a in 1:2^n_dim
      for b in 1:2^n_dim
        ∇ϕᵉ_a = H/detJ * ∇Φ(ξ, a)
        ∇ϕᵉ_b = H/detJ * ∇Φ(ξ, b)

        parcelaDerivada2 = α * dot(∇ϕᵉ_b, ∇ϕᵉ_a) * detJ
        @inbounds parcelaNormal = β * ϕᵉ[a] * ϕᵉ[b] * detJ

        @inbounds parcelaDerivada1 = 0 #γ * vec_Φ[a] * ∇ϕᵉ_b

        @inbounds Kᵉ[a,b] += WW * (parcelaDerivada2 + parcelaNormal + parcelaDerivada1)
      end
    end
  end
end

function montaK_geral(run_values::RunValues, malha::Malha)
  (; ne, neq, dx, EQ, LG, n_dim, coords, Nx, base) = malha
  X = coords

  npg = 2
  
  ϕξ, P, W = quadratura_ϕ(base, npg, n_dim)
  ∇ϕξ, P, W = quadratura_∇ϕ(base, npg, n_dim)
  
  Kᵉ = zeros(Float64, 2^n_dim, 2^n_dim)
  band = Nx[1]
  K = BandedMatrix(Zeros(neq, neq), (band, band))
  
  for e in 1:ne
    idx = LG[:,e]
    
    # Coordenadas dos vértices do elemento finito Ωᵉ
    Xe = ()
    for d in 1:n_dim
      Xe = (Xe..., X[d][idx])
    end
    
    idx = EQ[idx]

    montaKᵉ_geral!(Kᵉ, Xe, P, W, ϕξ, ∇ϕξ, n_dim, dx, run_values)

    for b in 1:2^n_dim
      @inbounds j = idx[b]
      for a in 1:2^n_dim
        @inbounds i = idx[a]
        if i <= neq && j <= neq
          @inbounds K[i,j] += Kᵉ[a,b]
        end
      end
    end
  end

  return sparse(K)
end

function montaFᵉ_geral!(Fᵉ, f, Xe, P, W, ϕξ, ∇ϕξ, n_dim)
  fill!(Fᵉ, 0.0)

  npg = length(P)
  for ξ in 1:npg
    ϕᵉ = ϕξ[1][ξ,:]

    x = ()
    for d in 1:n_dim
      x = (x..., dot(Xe[d], ϕᵉ))
    end
    
    M, detJ = jacobiano(n_dim, Xe, ∇ϕξ, ξ)
    @assert detJ > 0 "O determinante jacobiano deve ser positivo"

    WW = prod(W[ξ])
    for a in 1:2^n_dim
      @inbounds Fᵉ[a] += WW * f(x...) * ϕᵉ[a] * detJ
    end
  end
end

function montaF_geral(run_values::RunValues, malha::Malha)
  (;f) = run_values
  (; ne, neq, coords, LG, EQ, n_dim, base) = malha
  X = coords
  
  npg = 5

  ϕξ, P, W = quadratura_ϕ(base, npg, n_dim)
  ∇ϕξ, P, W = quadratura_∇ϕ(base, npg, n_dim)
  
  F = zeros(neq+1)
  Fᵉ = zeros(2^n_dim)

  for e in 1:ne
    idx = LG[:,e]
    
    # Coordenadas dos vértices do elemento finito Ωᵉ
    Xe = ()
    for d in 1:n_dim
      Xe = (Xe..., X[d][idx])
    end
    
    idx = EQ[idx]
    
    montaFᵉ_geral!(Fᵉ, f, Xe, P, W, ϕξ, ∇ϕξ, n_dim)
    
    for a in 1:2^n_dim
      F[idx[a]] += Fᵉ[a]
    end
  end
  
  xPTne = zeros(npg, ne)
  return F[1:neq], xPTne
end

function solveSys_geral(run_values::RunValues, malha::Malha)
  K = montaK_geral(run_values, malha)

  F, xPTne = montaF_geral(run_values, malha)

  C = zeros(Float64, malha.neq)

  C .= K\F

  F = nothing; K = nothing;

  return C, malha.EQoLG, xPTne
end