function avaliar_quadratura(base_func::Function, npg::Int64, n_funcs::Int64, n_dim::Int64)
  P, W = legendre(npg)

  phiP = zeros(npg, n_funcs)
  for a in 1:n_funcs
    for ksi in 1:npg
      phiP[ksi, a] = base_func(P[ksi])[a]
    end
  end
  return phiP, P, W
end

function xksi(ksi, e, X)
  h = X[e+1] - X[e]
  return h/2*(ksi+1) + X[e]
end

function montaK(run_values::RunValues, malha)
  (; alpha, beta, gamma) = run_values
  (; ne, neq, dx, EQoLG) = malha

  npg = 2
  
  phiP, P, W = avaliar_quadratura(ϕ_1D, npg, 2, 1)
  dphiP, P, W = avaliar_quadratura(∇ϕ_1D, npg, 2, 1)

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

function montaF(run_values::RunValues, malha::Malha)
  (; dx, ne, neq, coords, EQoLG) = malha
  (; f) = run_values
  
  X = coords[1]
  npg = 5

  phiP, P, W = avaliar_quadratura(ϕ_1D, npg, 2, 1)
  
  F = zeros(neq+1)
  xPTne = zeros(npg, ne)
  for e in 1:ne
    for ksi in 1:npg
      @inbounds fxptne = f(xksi(P[ksi], e, X))
      @inbounds xPTne[ksi, e] = fxptne
      for a in 1:2
        @inbounds i = EQoLG[a, e]

        @inbounds F[i] += W[ksi] * fxptne * phiP[ksi, a] * dx/2
      end
    end
  end

  return F[1:neq], xPTne
end

function solveSys(run_values::RunValues, malha::Malha)
  K = montaK(run_values, malha)

  F, xPTne = montaF(run_values, malha)

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