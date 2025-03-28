function montaLG(ne, base)
  LG = zeros(Int64, ne, base.nB)

  LG1 = [(base.nB-1)*i - (base.nB-2) for i in 1:ne]
  
  for e in 1:ne
      for a in 1:base.nB
          @inbounds LG[e,a] = LG1[e] + (a-1)
      end
  end

  return LG'
end

function montaEQ(ne, neq, base)
  EQ = zeros(Int64, (base.nB-1)*ne+1, 1)

  EQ[1] = neq+1; EQ[end] = neq+1
  
  for i in 2:(base.nB-1)*ne
      @inbounds EQ[i] = i-1
  end

  return EQ
end

function PHI(P, n_dim=1)
  phis = 2.0^(-n_dim) * ones(Float64, 2^n_dim)
  for d in 1:n_dim
    for n_phi in 1:2^n_dim
        sinal = floor((2.0^(d-1) + n_phi-1) * 2.0^-(d-1))
        phis[n_phi] *= 1 + ((-1.0)^sinal)*P[d]
    end
  end

  return phis
end

function dPHI(P, n_dim=1)
  if n_dim == 1
    return 2.0^-n_dim * [
      1-P[1],
      1+P[1]
    ]
  elseif n_dim == 2
    return 2.0^-n_dim * [
      (1-P[1]) * (1-P[2]),
      (1+P[1]) * (1-P[2]),
      (1-P[1]) * (1+P[2]),
      (1+P[1]) * (1+P[2]),
    ]
  else
    return Nothing
  end
end

function xksi(ksi, e, X)
  h = X[e+1] - X[e]
  return h/2*(ksi+1) + X[e]
end


function montaK(base, ne, neq, dx, alpha, beta, gamma, sigma, EQoLG::Matrix{Int64}, X)
  npg = base.p+1; P, W = legendre(npg)

  phiP(ksi, a) = PHI(ksi, length(ksi))[a]; dphiP(ksi, a) = dPHI(ksi, length(ksi))[a];

  Ke = zeros(Float64, 2, 2)

  for a in 1:2
      for b in 1:2
          for ksi in 1:npg
              @inbounds parcelaNormal = beta*dx/2 * W[ksi] * phiP(P[ksi], a) * phiP(P[ksi], b);
              @inbounds parcelaDerivada1 = gamma * W[ksi] * phiP(P[ksi], a) * dphiP(P[ksi], b);
              @inbounds parcelaDerivada2 = 2*alpha/dx * W[ksi] * dphiP(P[ksi], a) * dphiP(P[ksi], b);

              @inbounds Ke[a,b] += parcelaDerivada2 + parcelaNormal + parcelaDerivada1
          end
      end
  end

  K = BandedMatrix(Zeros(neq, neq), (base.nB-1, base.nB-1))
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

function montaF(base, ne, neq, X, f::Function, EQoLG)
  npg = 5; P, W = legendre(npg)

  phiP = zeros(npg, base.nB)
  for a in 1:2
      for ksi in 1:npg
          phiP[ksi, a] = PHI(P[ksi], length(P[ksi]))[a]
      end
  end
  
  dx = 1/ne
  
  F = zeros(neq+1)
  xPTne = zeros(npg, ne)
  for e in 1:ne
      for ksi in 1:npg
          @inbounds fxptne = f(xksi(P[ksi], e, X))
          @inbounds xPTne[ksi, e] = fxptne
          for a in 1:2
              @inbounds partial = dx/2 * W[ksi] * phiP[ksi, a]
              @inbounds i = EQoLG[a, e]
              @inbounds F[i] += partial * fxptne
          end
      end
  end

  return F[1:neq], xPTne
end

function solveSys(base, alpha, beta, gamma, sigma, ne, a, b, f, u)
  dx = 1/ne;
  neq = base.neq;

  X = a:dx:b
  
  EQ = montaEQ(ne, neq, base); LG = montaLG(ne, base)
  EQoLG = EQ[LG]

  EQ = nothing; LG = nothing;
  
  K = montaK(base, ne, neq, dx, alpha, beta, gamma, sigma, EQoLG, X)

  F, xPTne = montaF(base, ne, neq, X, f, EQoLG)

  C = zeros(Float64, neq)

  C .= K\F

  F = nothing; K = nothing;

  return C, EQoLG, xPTne
end