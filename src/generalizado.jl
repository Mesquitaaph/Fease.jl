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

function avaliar_quadratura_geral_ϕ(base_func::Function, npg::Int64, n_funcs::Int64, n_dim::Int64)
  P, W = quadratura_gauss(npg, n_dim)

  ϕPₙ = ()
  for d in 1:n_dim
    ϕP = zeros(npg^n_dim, n_funcs^n_dim)
    for ξ in 1:npg^n_dim
      ϕP[ξ, :] .= base_func(P[ξ]...)[1]
    end
    ϕPₙ = (ϕPₙ..., ϕP)
  end

  return ϕPₙ, P, W
end

function avaliar_quadratura_geral(base_func::Function, npg::Int64, n_funcs::Int64, n_dim::Int64)
  P, W = quadratura_gauss(npg, n_dim)

  ϕPₙ = ()
  for d in 1:n_dim
    ϕP = zeros(npg^n_dim, n_funcs^n_dim)
    for ξ in 1:npg^n_dim
      ϕP[ξ, :] .= base_func(P[ξ]...)[d]
    end
    ϕPₙ = (ϕPₙ..., ϕP)
  end

  return ϕPₙ, P, W
end

function montaK_geral(run_values::RunValues, malha)
  (; alpha, beta, gamma) = run_values
  (; ne, neq, dx, EQoLG) = malha

  npg = 2
  
  phiP, P, W = avaliar_quadratura_geral(ϕ_geral, npg, 2, 1)
  dphiP, P, W = avaliar_quadratura_geral(∇ϕ_1D, npg, 2, 1)

  Ke = zeros(Float64, 2, 2)
  for a in 1:2
      for b in 1:2
          for ksi in 1:npg
              WW = prod(W[ksi])
              @inbounds parcelaNormal = beta*dx/2 * WW * phiP[ksi, a] * phiP[ksi, b];
              @inbounds parcelaDerivada1 = gamma * WW * phiP[ksi, a] * dphiP[ksi, b];
              @inbounds parcelaDerivada2 = 2*alpha/dx * WW * dphiP[ksi, a] * dphiP[ksi, b];

              @inbounds Ke[a,b] += parcelaDerivada2 + parcelaNormal + parcelaDerivada1
          end
      end
  end

  K = BandedMatrix(Zeros(neq, neq), (1, 1))
  for e in 1:ne
      for b in 1:2
          @inbounds j = EQoLG[b, e]
          for a in 1:2
              @inbounds i = EQoLG[a, e]
              if i <= neq && j <= neq
                  @inbounds K[i,j] += Ke[a,b]
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
    vec_∂ϕ_∂ξ₁ = ∇ϕξ[1][ξ,:]

    vec_ϕ = ϕξ[1][ξ,:]

    x = ()
    for d in 1:n_dim
      x = (x..., dot(Xe[d], vec_ϕ))
    end
    
    M_J = Matrix{Float64}(undef, n_dim, n_dim)
    for i in 1:n_dim
      for j in 1:n_dim
        M_J[i,j] = dot(Xe[i], ∇ϕξ[j][ξ,:])
      end
    end

    detJ = det(M_J)
    @assert detJ > 0 "O determinante jacobiano deve ser positivo"

    WW = prod(W[ξ])
    for a in 1:2^n_dim
      @inbounds Fᵉ[a] += WW * f(x...) * vec_ϕ[a] * detJ
    end
  end
end

function montaF_geral(run_values::RunValues, malha::Malha)
  (;f) = run_values
  (; ne, neq, coords, LG, EQ, n_dim) = malha
  X = coords
  
  npg = 5

  ϕξ, P, W = avaliar_quadratura_geral_ϕ(ϕ_geral, npg, 2, n_dim)
  ∇ϕξ, P, W = avaliar_quadratura_geral(∇ϕ_geral, npg, 2, n_dim)
  
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