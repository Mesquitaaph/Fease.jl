function avaliar_quadratura_geral(base_func::Function, npg::Int64, n_funcs::Int64, n_dim::Int64)
  P, W = legendre(npg)

  phiP = zeros(npg, n_funcs)
  for a in 1:n_funcs
      for ksi in 1:npg
          phiP[ksi, a] = base_func(P[ksi], n_dim)[a]
      end
  end
  return phiP, P, W
end

function xksi_geral(ksi, e, X)
  h = X[e+1] - X[e]
  return h/2*(ksi+1) + X[e]
end

function montaK_geral(run_values::RunValues, malha)
  (; alpha, beta, gamma) = run_values
  (; ne, neq, dx, EQoLG) = malha

  npg = 2
  
  phiP, P, W = avaliar_quadratura_geral(ϕ_1D, npg, 2, 1)
  dphiP, P, W = avaliar_quadratura_geral(∇ϕ_1D, npg, 2, 1)

  Ke = zeros(Float64, 2, 2)
  for a in 1:2
      for b in 1:2
          for ksi in 1:npg
              @inbounds parcelaNormal = beta*dx/2 * W[ksi] * phiP[ksi, a] * phiP[ksi, b];
              @inbounds parcelaDerivada1 = gamma * W[ksi] * phiP[ksi, a] * dphiP[ksi, b];
              @inbounds parcelaDerivada2 = 2*alpha/dx * W[ksi] * dphiP[ksi, a] * dphiP[ksi, b];

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

function montaFᵉ_geral!(Fᵉ, f, X1e, P, W, ϕξ, ∇ϕξ)
  fill!(Fᵉ, 0.0)

  npg = length(P)
  for i in 1:npg
    vec_ϕ = ϕξ[i,:]

    x = dot(X1e, vec_ϕ)

    vec_∂ϕ_∂ξ₁ = ∇ϕξ[i,:]
    
    detJ = dot(X1e, vec_∂ϕ_∂ξ₁)
    @assert detJ > 0 "O determinante jacobiano deve ser positivo"

    for a in 1:2
      @inbounds Fᵉ[a] += W[i] * f(x) * ϕξ[i, a] * detJ
    end
  end
end

function montaF_geral(f, malha::Malha)
  (; ne, neq, coords, LG, EQ, EQoLG) = malha
  X = coords[1]
  
  npg = 5

  ϕξ, P, W = avaliar_quadratura_geral(ϕ_1D, npg, 2, 1)
  ∇ϕξ, P, W = avaliar_quadratura_geral(∇ϕ_1D, npg, 2, 1)
  
  F = zeros(neq+1)
  Fᵉ = zeros(4)

  for e in 1:ne
    idx = LG[:,e]
    
    # Coordenadas dos vértices do elemento finito Ωᵉ
    X1e = X[idx]
    # X2e = X₂[idx]
    
    idx = EQ[idx]
    
    montaFᵉ_geral!(Fᵉ, f, X1e, P, W, ϕξ, ∇ϕξ)
    
    for a in 1:2
      F[idx[a]] += Fᵉ[a]
    end
  end
  
  xPTne = zeros(npg, ne)
  return F[1:neq], xPTne
end

function solveSys_geral(run_values::RunValues, malha::Malha)
  K = montaK_geral(run_values, malha)

  F, xPTne = montaF_geral(run_values.f, malha)

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