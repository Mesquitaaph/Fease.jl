using Fease

function malhas()
  # Define o número de sub-intervalos no eixo horizontal e vertical, respectivamente.
  Nx1, Nx2 = 4, 4

  # Define o tipo da base de funções interpoladoras do subespaço aproximado Vₘ.
  baseType = BaseTypes.linearLagrange

  # Define extremos do intervalo da malha.
  a = (0.0, 0.0)
  b = (1.0, 1.0)

  # Define e plota a malha uniforme com os valores atribuídos acima.
  malha = monta_malha_2D_uniforme(baseType, Nx1, Nx2, a, b)
  plot_malha_2D(malha)

  # Define e plota a malha uniforme, com ruído, utilizando os valores atribuídos acima.
  malha = malha2D_adiciona_ruido(malha)
  plot_malha_2D(malha)

  # Ponto foco
  # Precisao 0
  ponto_foco = (0.375, 0.375)
  precisao = 0

  # Função que constrói uma malha 2D que cria mais elementos próximos ao ponto_foco, com 0 de precisão.
  malha = monta_malha_2D_foco(baseType, Nx1, Nx2, a, b, ponto_foco, precisao)
  plot_malha_2D(malha)

  # Precisao 1
  precisao = 1

  # Função que constrói uma malha 2D que cria mais elementos próximos ao ponto_foco, com 1 de precisão.
  malha = monta_malha_2D_foco(baseType, Nx1, Nx2, a, b, ponto_foco, precisao)
  plot_malha_2D(malha)

  # Precisao 2
  precisao = 2

  # Função que constrói uma malha 2D que cria mais elementos próximos ao ponto_foco, com 2 de precisão.
  malha = monta_malha_2D_foco(baseType, Nx1, Nx2, a, b, ponto_foco, precisao)
  return plot_malha_2D(malha)
end

malhas()
