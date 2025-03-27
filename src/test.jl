function single_run()
  println("Rodando \n")
  @benchmark begin
      alpha, beta, gamma, sigma, a, b, u, u_x, f = examples(3); ne = 2^2; baseType = BaseTypes.linearLagrange
      base = LocalBases[Symbol(baseType)](ne)
      C, EQoLG, xPTne = solveSys(base, alpha, beta, gamma, sigma, ne, a, b, f, u)

  end
end

P2 = [3, 7, 23]
n_dim = length(P2)

b2_1 = @benchmark PHI(P2, n_dim)
b2_2 = @benchmark PHI_original(P2, n_dim)

PHI(P2, n_dim)