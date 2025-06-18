using MyProject

α, β, f, u = exemplo1()

run_values = RunValues(α, β, 0.0, f, u)

Nx1, Nx2 = 4, 3;
ne = Nx1 * Nx2

baseType = BaseTypes.linearLagrange
base = monta_base(baseType, ne)

a = (0.0, 0.0);
b = (1.0, 1.0)
malha = monta_malha_2D_uniforme(base, Nx1, Nx2, a, b)

K = montaK_geral(run_values, malha)

F = montaF_geral(run_values, malha)

C = K\F

display(C)
