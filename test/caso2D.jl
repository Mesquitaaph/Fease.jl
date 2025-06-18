function single_run_2D(ruido::Bool = false)
  α, β, f, u = exemplo1()

  run_values = RunValues(α, β, 0.0, f, u)

  Nx1, Nx2 = 4, 3;
  ne = Nx1 * Nx2

  baseType = BaseTypes.linearLagrange
  base = monta_base(baseType, ne)

  a = (0.0, 0.0);
  b = (1.0, 1.0)
  malha = monta_malha_2D_uniforme(base, Nx1, Nx2, a, b)

  if ruido
    # Exemplo 2: Adiciona ruído nos nós internos
    Random.seed!(42) # Define uma semente para reprodutibilidade
    malha = malha2D_adiciona_ruido(malha)
  else
    # Exemplo 1: Malha uniforme de retângulos
  end

  K = montaK_geral(run_values, malha)

  F = montaF_geral(run_values, malha)

  C = K\F

  return C
end

@testset "caso2D.jl" begin
  correct_result_ex1 = [0.659197679603509
                        0.9322462987801567
                        0.659197679603509
                        0.6591976796035091
                        0.9322462987801565
                        0.659197679603509]
  ruido = false
  @test single_run_2D(ruido) ≈ correct_result_ex1

  correct_result_ex2 = [0.7065102277637031
                        0.969014464383587
                        0.6584488958892918
                        0.7346980208325475
                        0.8893466262809584
                        0.7082557097977525]
  ruido = true
  @test single_run_2D(ruido) ≈ correct_result_ex2
end
