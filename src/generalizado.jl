function avaliar_quadratura_geral(base_func::Function, npg::Int64, n_funcs::Int64, n_dim::Int64)
  p, w = legendre(npg)

  P = Vector{Tuple}(undef, npg^n_dim)
  W = similar(P)

  for i in 1:npg^n_dim
    ponto = ()
    peso = ()
    for d in 1:n_dim
      idx_global = floor((2.0^(n_dim-d) + i-1) * 2.0^-(n_dim-d))
      idx = Int((idx_global-1) % npg)+1
      ponto = (ponto..., p[idx])
      peso = (peso..., w[idx])
    end
    P[i] = ponto
    W[i] = peso
  end

  ϕP = zeros(npg^n_dim, n_funcs^n_dim)
  for ξ in 1:npg^n_dim
    ϕP[ξ, :] .= base_func(P[ξ]...)
  end

  return ϕP, P, W
end

function xksi_geral(ksi, e, X)
  h = X[e+1] - X[e]
  return h/2*(ksi+1) + X[e]
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
  X1e = Xe[1]

  npg = length(P)
  for i in 1:npg
    vec_∂ϕ_∂ξ₁ = ∇ϕξ[i,:]

    vec_ϕ = ϕξ[i,:]

    x = dot(X1e, vec_ϕ)
    
    detJ = dot(X1e, vec_∂ϕ_∂ξ₁)
    @assert detJ > 0 "O determinante jacobiano deve ser positivo"

    WW = prod(W[i])
    for a in 1:2^n_dim
      @inbounds Fᵉ[a] += WW * f(x) * ϕξ[i, a] * detJ
    end
  end
end

function montaF_geral(run_values::RunValues, malha::Malha)
  (;f) = run_values
  (; ne, neq, coords, LG, EQ, n_dim) = malha
  X = coords
  
  npg = 5

  ϕξ, P, W = avaliar_quadratura_geral(ϕ_geral, npg, 2, 1)
  ∇ϕξ, P, W = avaliar_quadratura_geral(∇ϕ_1D, npg, 2, 1)
  
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

  # display("Matriz K RAPHAEL:")
  # display(K)
  # display("Vetor F RAPHAEL:")
  # display(F)

  C = zeros(Float64, malha.neq)

  C .= K\F
  # display("Solução aproximada U RAPHAEL:")
  # display(C)

  F = nothing; K = nothing;

  return C, malha.EQoLG, xPTne
end