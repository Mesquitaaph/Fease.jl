function erro_L2(base, ne, EQoLG, C, u, u_x, xPTne)
  npg = 5; P, W = legendre(npg)
  
  phiP(ksi, a) = PHI(ksi, base)[a]; dphiP(ksi, a) = dPHI(ksi, base)[a];

  h = 1/ne

  cEQoLG = vcat(C, 0)[EQoLG]

  EL2 = 0.0
  for e in 1:ne
      @inbounds c1e = cEQoLG[1, e]
      @inbounds c2e = cEQoLG[2, e]
      for ksi in 1:npg
         @inbounds EL2 += W[ksi] * (u(xPTne[ksi, e]) - (phiP(P[ksi], 1) * c1e) - (phiP(P[ksi], 2) * c2e))^2
      end
  end

  EL2 = sqrt(h/2 * EL2)

  EH01 = 0.0
  # EH01 = sqrt(h/2 * sum(W' * ((u_x.(xPTne) - (2/h .* dphiP * cEQoLG)).^2)))

  return EL2, EH01
end

function convergence_test!(NE, E, dE, example)
  alpha, beta, gamma, sigma, a, b, u, u_x, f = examples(example)

  for i = 1:lastindex(NE)
      base = LocalBases[Symbol(baseType)](NE[i])
      Ci, EQoLGi, xPTnei = solveSys(base, alpha, beta, gamma, sigma, NE[i], a, b, f, u)
      E[i], dE[i] = erro_L2(base, NE[i], EQoLGi, Ci, u, u_x, xPTnei)
  end
end