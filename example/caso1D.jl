using MyProject

example = 1

alpha, beta, gamma, a, b, u, u_x, f = examples(example); Nx1 = ne = 2^3

run_values = RunValues(alpha, beta, gamma, f, u)

baseType = BaseTypes.linearLagrange
base = monta_base(baseType, ne)

malha = monta_malha_1D_uniforme(Nx1, base, a, b)

K = montaK_geral(run_values, malha)

F = montaF_geral(run_values, malha)

C = K\F

display(C)