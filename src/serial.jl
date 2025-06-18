function avaliar_quadratura(base_func::Function, npg::Int64, n_funcs::Int64, n_dim::Int64)
  P, W = legendre(npg)

  phiP = zeros(npg, n_funcs)
  for a in 1:n_funcs
    for ksi in 1:npg
      phiP[ksi, a] = base_func(P[ksi])[1][a]
    end
  end
  return phiP, P, W
end

function xksi(ksi, e, X)
  h = X[e + 1] - X[e]
  return h/2*(ksi+1) + X[e]
end

function montaK_1D(run_values::RunValues, malha)
  (; α, β, γ) = run_values
  (; ne, neq, dx, EQoLG) = malha

  npg = 2

  phiP, P, W = avaliar_quadratura(ϕ_1D, npg, 2, 1)
  dphiP, P, W = avaliar_quadratura(∇ϕ_1D, npg, 2, 1)

  Ke = zeros(Float64, 2, 2)
  @inbounds for a in 1:2
    for b in 1:2
      for ksi in 1:npg
        parcelaNormal = β*dx/2 * W[ksi] * phiP[ksi, a] * phiP[ksi, b];
        parcelaDerivada1 = γ * W[ksi] * phiP[ksi, a] * dphiP[ksi, b];
        parcelaDerivada2 = 2*α/dx * W[ksi] * dphiP[ksi, a] * dphiP[ksi, b];

        Ke[a, b] += parcelaDerivada2 + parcelaNormal + parcelaDerivada1
      end
    end
  end

  K = BandedMatrix(Zeros(neq, neq), (1, 1))
  @inbounds for e in 1:ne
    for b in 1:2
      j = EQoLG[b, e]
      for a in 1:2
        i = EQoLG[a, e]
        if i <= neq && j <= neq
          K[i, j] += Ke[a, b]
        end
      end
    end
  end

  return sparse(K)
end

function montaF_1D(run_values::RunValues, malha::Malha)
  (; dx, ne, neq, coords, EQoLG) = malha
  (; f) = run_values

  X = coords[1]
  npg = 5

  phiP, P, W = avaliar_quadratura(ϕ_1D, npg, 2, 1)

  F = zeros(neq+1)
  xPTne = zeros(npg, ne)
  @inbounds for e in 1:ne
    for ksi in 1:npg
      fxptne = f(xksi(P[ksi], e, X))
      xPTne[ksi, e] = fxptne
      for a in 1:2
        i = EQoLG[a, e]

        F[i] += W[ksi] * fxptne * phiP[ksi, a] * dx/2
      end
    end
  end

  return F[1:neq], xPTne
end

function solveSys_1D(run_values::RunValues, malha::Malha)
  K = montaK_1D(run_values, malha)

  F, xPTne = montaF_1D(run_values, malha)

  C = zeros(Float64, malha.neq)

  C .= K\F

  F = nothing;
  K = nothing;

  return C, malha.EQoLG, xPTne
end
