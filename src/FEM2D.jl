function x₁_de_ξ(ξ₁::Float64, h₁::Float64, p₁::Float64)::Float64
  x₁ = p₁ + (h₁ / 2) * (ξ₁ + 1)
  return x₁
end

function x₂_de_ξ(ξ₂::Float64, h₂::Float64, p₂::Float64)::Float64
  x₂ = p₂ + (h₂ / 2) * (ξ₂ + 1)
  return x₂
end

############### MONTA K ###############
function monta_Kᵉ_quadrilatero!(
    Kᵉ::Matrix{Float64}, run_values::RunValues, X1e::Vector{Float64},
    X2e::Vector{Float64}, P::Vector{Float64}, W::Vector{Float64})
  # Zera as entradas da matriz local Kᵉ
  fill!(Kᵉ, 0.0)

  (; α, β) = run_values

  # Loop de quadratura dupla (ξ₁, ξ₂) para integração gaussiana
  for i in 1:length(P)  # Pontos de quadratura no eixo ξ₁
    ξ₁ = P[i]

    # Calcula as derivadas de ϕ em relação a ξ₂ no ponto (ξ₁,ξ₂)
    # Obs.: As funções em ∂ϕ_∂ξ₂ dependem apenas de ξ₁
    vec_∂ϕ_∂ξ₂ = ∂ϕ_∂ξ₂(ξ₁)

    for j in 1:length(P)  # Pontos de quadratura no eixo ξ₂
      ξ₂ = P[j]

      # Vetor com o valor em (ξ₁,ξ₂) das funções ϕ₁, ϕ₂, ϕ₃ e ϕ₄
      vec_ϕ = ϕ_2D(ξ₁, ξ₂)

      # Calcula as derivadas de ϕ em relação a ξ₁ no ponto (ξ₁,ξ₂)
      # Obs.: As funções em ∂ϕ_∂ξ₁ dependem apenas de ξ₂
      vec_∂ϕ_∂ξ₁ = ∂ϕ_∂ξ₁(ξ₂)

      # Calcula as entradas da matriz Jacobiana
      ∂x₁_∂ξ₁ = dot(X1e, vec_∂ϕ_∂ξ₁)
      ∂x₁_∂ξ₂ = dot(X1e, vec_∂ϕ_∂ξ₂)
      ∂x₂_∂ξ₁ = dot(X2e, vec_∂ϕ_∂ξ₁)
      ∂x₂_∂ξ₂ = dot(X2e, vec_∂ϕ_∂ξ₂)

      # Calcula o determinante Jacobiano do mapeamento isoparamétrico
      detJ = ∂x₁_∂ξ₁ * ∂x₂_∂ξ₂ - ∂x₁_∂ξ₂ * ∂x₂_∂ξ₁
      @assert detJ>0 "O determinante jacobiano deve ser positivo"

      # Calcula as entradas da matriz Hᵀ * H
      HᵀH₁₁ = ∂x₂_∂ξ₂^2 + ∂x₁_∂ξ₂^2
      HᵀH₁₂ = -∂x₁_∂ξ₁ * ∂x₁_∂ξ₂ - ∂x₂_∂ξ₁ * ∂x₂_∂ξ₂
      HᵀH₂₂ = ∂x₂_∂ξ₁^2 + ∂x₁_∂ξ₁^2

      # Calcula a contribuição de quadratura e acumula o valor na matriz Kᵉ
      for b in 1:4
        for a in 1:4
          Kᵉ[a,
          b] += W[i] * W[j] *
                (
                  (α / detJ) *
                  (
                  vec_∂ϕ_∂ξ₁[b] * (HᵀH₁₁ * vec_∂ϕ_∂ξ₁[a] + HᵀH₁₂ * vec_∂ϕ_∂ξ₂[a]) +
                  vec_∂ϕ_∂ξ₂[b] * (HᵀH₁₂ * vec_∂ϕ_∂ξ₁[a] + HᵀH₂₂ * vec_∂ϕ_∂ξ₂[a])
                ) + β * vec_ϕ[b] * vec_ϕ[a] * detJ
                )
        end
      end
    end
  end
end

function monta_K_quadrilatero(
    run_values::RunValues, malha::Malha)::SparseMatrixCSC{Float64, Int64}
  (; coords, neq, EQ, LG) = malha
  X₁, X₂ = coords
  m = neq

  # Número de elementos na malha
  ne = size(LG, 2)

  # Pontos e pesos de quadratura de Gauss-Legendre
  P, W = legendre(2)

  # Inicializa a matriz local Kᵉ
  Kᵉ = zeros(4, 4)

  # Inicializa a matriz esparsa K com tamanho (m+1) x (m+1)
  K = spzeros(m + 1, m + 1)

  # Loop sobre os elementos
  for e in 1:ne
    # Índices dos vértices/funções globais do elemento finito Ωᵉ
    idx = LG[:, e]

    # Coordenadas dos vértices do elemento finito Ωᵉ
    X1e = X₁[idx]
    X2e = X₂[idx]

    # Numeração de equação dos índices no vetor idx
    idx = EQ[idx]

    # Calcula a matriz local Kᵉ
    monta_Kᵉ_quadrilatero!(Kᵉ, run_values, X1e, X2e, P, W)

    # Loop sobre as colunas (b) e linhas (a) da matriz local Kᵉ
    for b in 1:4
      j = idx[b]
      for a in 1:4
        i = idx[a]
        K[i, j] += Kᵉ[a, b]
      end
    end
  end

  # Remove a última linha e coluna da matriz K
  return K[1:m, 1:m]
end

############### MONTA F ###############
function monta_Fᵉ_quadrilatero!(Fᵉ::Vector{Float64}, f::Function, X1e::Vector{Float64},
    X2e::Vector{Float64}, P::Vector{Float64}, W::Vector{Float64})
  # Zera as entradas do vetor local Fᵉ
  fill!(Fᵉ, 0.0)

  # Integração via quadratura gaussiana dupla
  for i in 1:length(P)        # Loop sobre os pontos de quadratura no eixo ξ₁
    ξ₁ = P[i]

    # Calcula as derivadas de ϕ em relação a ξ₂ no ponto (ξ₁,ξ₂)
    # Obs.: As funções em ∂ϕ_∂ξ₂ dependem apenas de ξ₁
    vec_∂ϕ_∂ξ₂ = ∂ϕ_∂ξ₂(ξ₁)

    for j in 1:length(P)      # Loop sobre os pontos de quadratura no eixo ξ₂
      ξ₂ = P[j]

      # Vetor com o valor em (ξ₁,ξ₂) das funções ϕ₁, ϕ₂, ϕ₃ e ϕ₄
      vec_ϕ = ϕ_2D(ξ₁, ξ₂)

      # Calcula as derivadas de ϕ em relação a ξ₁ no ponto (ξ₁,ξ₂)
      # Obs.: As funções em ∂ϕ_∂ξ₁ dependem apenas de ξ₂
      vec_∂ϕ_∂ξ₁ = ∂ϕ_∂ξ₁(ξ₂)

      # Coordenadas físicas (x₁, x₂) mapeadas a partir de (ξ₁, ξ₂)
      x₁ = dot(X1e, vec_ϕ)
      x₂ = dot(X2e, vec_ϕ)

      # Calcula o determinante jacobiano do mapeamento isoparamétrico
      detJ = dot(X1e, vec_∂ϕ_∂ξ₁) * dot(X2e, vec_∂ϕ_∂ξ₂) -
             dot(X1e, vec_∂ϕ_∂ξ₂) * dot(X2e, vec_∂ϕ_∂ξ₁)
      @assert detJ>0 "O determinante jacobiano deve ser positivo"

      # Calcula a contribuição de quadratura e acumula em Fᵉ
      for a in 1:4
        Fᵉ[a] += W[i] * W[j] * f(x₁, x₂) * vec_ϕ[a] * detJ
      end
    end
  end
end

function monta_F_quadrilatero(run_values::RunValues, malha::Malha)::Vector{Float64}
  (; f) = run_values

  X₁, X₂ = malha.coords

  m, EQ = malha.neq, malha.EQ
  LG = malha.LG

  # Número total de elementos finitos na malha
  ne = malha.ne

  # Pontos e pesos de quadratura de Gauss-Legendre de ordem 5
  P, W = legendre(5)

  # Inicializa os vetores locais e globais
  Fᵉ = zeros(4)
  F = zeros(m + 1)

  # Itera sobre os elementos Ωᵉ
  for e in 1:ne
    # Índices dos vértices/funções globais do elemento finito Ωᵉ
    idx = LG[:, e]

    # Coordenadas dos vértices do elemento finito Ωᵉ
    X1e = X₁[idx]
    X2e = X₂[idx]

    # Numeração de equação dos índices no vetor idx
    idx = EQ[idx]

    # Monta o vetor local Fᵉ
    monta_Fᵉ_quadrilatero!(Fᵉ, f, X1e, X2e, P, W)

    # Adiciona a contribuição do elemento finito ao vetor global F
    for a in 1:4
      F[idx[a]] += Fᵉ[a]
    end
  end

  # Retorna o vetor global F com tamanho `m`, excluindo a última entrada adicional
  return F[1:m]
end

############### TESTES ###############
function plot_solução_aproximada(
    c̄::Vector{Float64}, Nx1::Int64, Nx2::Int64, EQoLG::Matrix{Int64})
  # Comprimentos da base (h₁) e altura (h₂) de cada elemento retangular Ωᵉ
  h₁, h₂ = 1 / Nx1, 1 / Nx2

  # Define uma discretização do intervalo de referência [-1, 1] nos eixos ξ₁ e ξ₂
  ξ₁ = collect(-1:0.1:1)
  ξ₂ = collect(-1:0.1:1)

  # Inicializa um objeto de gráfico que acumulará as superfícies
  plt = Plots.plot(seriestype = :surface, title = "Solução Aproximada",
    color = :viridis,
    xlabel = "x₁", ylabel = "x₂", colorbar = false, zlims = (0, 1),
    xticks = [0, 0.5, 1], yticks = [0, 0.5, 1], zticks = [0, 1])

  # Loop nos elementos da malha (percorrendo cada subdivisão ao longo de x₂ e x₁)
  for j in 1:Nx2
    # Define a segunda coordenada do ponto inferior esquerdo do retângulo `Ωᵉ`.
    p₂ = (j - 1) * h₂
    # Calcula os valores correspondentes no eixo x₂ para os pontos ξ₂ em [-1,1]
    x₂ = x₂_de_ξ.(ξ₂, h₂, p₂)

    for i in 1:Nx1
      # Primeira coordenada do ponto inferior esquerdo do retângulo `Ωᵉ`.
      p₁ = (i - 1) * h₁
      # Valores correspondentes no eixo x₁ para os pontos ξ₁ em [-1,1]
      x₁ = x₁_de_ξ.(ξ₁, h₁, p₁)

      # Determina a numeração do elemento finito atual (`e`) na malha
      e = (j - 1) * Nx1 + i

      # Obtém os coeficientes `c` da solução aproximada no elemento `e`
      ce = c̄[EQoLG[:, e]]

      # Solução aproximada na malha x₁ × x₂. Tamanho length(x₂) x length(x₁)
      mat = [dot(ce, ϕ_2D(ξ₁[a], ξ₂[b])) for b in 1:length(x₂), a in 1:length(x₁)]

      # Exibindo os resultados
      # display("------------------------------")
      # display("e = $e")
      # display("p1, p2 = $p₁, $p₂")
      # display("ce = $ce")
      # display("x₁ = $x₁")
      # display("x₂ = $x₂")
      # display("uₕ(x₁[a],x₂[b]) for b in 1:length(x₂), a in 1:length(x₁)")
      # display(mat)

      # Gera o gráfico da solução aproximada sobre o elemento `e`
      Plots.plot!(plt, x₁, x₂, mat, seriestype = :surface, alpha = 0.5)
    end
  end

  return plt
end

############### MALHA 2D ###############
function malha2D(Nx1, Nx2)
  # Define o comprimento da base (h₁) e altura (h₂) de cada elemento retangular Ωᵉ
  h₁, h₂ = 1 / Nx1, 0.5 / Nx2

  # Define a discretização em x₁ e x₂
  x₁ = collect(0:h₁:1)
  x₂ = collect(0:h₂:1)

  # Define as coordenadas de cada nó da malha
  X₁ = [x₁[i] for i in 1:(Nx1+1), j in 1:(Nx2+1)]
  X₂ = [x₂[j] for i in 1:(Nx1+1), j in 1:(Nx2+1)]
  # display(X₁)
  # display(X₂)

  return X₁, X₂, h₁, h₂
end

function malha2D_adiciona_ruido!(
    X₁::Matrix{Float64}, X₂::Matrix{Float64}, h₁::Float64, h₂::Float64)
  # Verifica se as matrizes têm a mesma dimensão
  @assert size(X₁)==size(X₂) "X₁ e X₂ devem ter as mesmas dimensões"

  # Verifica se a malha é suficientemente grande para ter nós internos
  if size(X₁, 1) > 2 && size(X₁, 2) > 2
    # Define os limites do ruído
    ruído_limite₁, ruído_limite₂ = h₁ / 4, h₂ / 4

    # Aplica ruído uniforme aos nós internos
    X₁[2:(end-1),
    2:(end-1)] .+= ruído_limite₁ *
                   (rand(Float64, size(X₁[2:(end-1), 2:(end-1)])) .- 0.5) * 2
    X₂[2:(end-1),
    2:(end-1)] .+= ruído_limite₂ *
                   (rand(Float64, size(X₂[2:(end-1), 2:(end-1)])) .- 0.5) * 2
  end
end

function plot_malha2D(
    X₁::AbstractArray{Float64}, X₂::AbstractArray{Float64}, LG::Matrix{Int64})
  # Cria uma nova figura
  fig = Plots.plot(legend = false, aspect_ratio = :equal,
    xticks = 0:0.25:1,  # Defina os ticks do eixo x
    yticks = 0:0.25:1)  # Defina os ticks do eixo y

  # Adiciona os nós da malha
  Plots.scatter!(X₁, X₂, markersize = 4, color = :blue)

  # Adiciona as linhas de contorno dos elementos
  for e in 1:size(LG, 2)
    # Obtém os índices dos nós do elemento `e`
    i₁, i₂, i₃, i₄ = LG[:, e]

    # Adiciona as linhas de contorno do elemento `e`
    Plots.plot!([X₁[i₁], X₁[i₂], X₁[i₄], X₁[i₃], X₁[i₁]],
      [X₂[i₁], X₂[i₂], X₂[i₄], X₂[i₃], X₂[i₁]], color = :black)
  end

  return fig
end

function exemplo_malha2D(Nx1, Nx2; ruido = true)
  # Cria a malha original
  X₁, X₂, h₁, h₂ = malha2D(Nx1, Nx2)

  # Adiciona ruido
  if ruido
    malha2D_adiciona_ruido!(X₁, X₂, h₁, h₂)
  end
  # Define a matriz de conectividade local/global (LG)
  LG = monta_LG_2D(Nx1, Nx2)

  # Cria a figura da malha 2D
  fig = plot_malha2D(X₁, X₂, LG)

  return fig
end

let
  Nx1 = 4
  Nx2 = 3

  fig = exemplo_malha2D(Nx1, Nx2; ruido = false)  # Aqui a função é chamada inicialmente
  xlabel!(fig, "x_1")
  ylabel!(fig, "x_2")
  yticks!(fig, 0:(1/Nx2):1)             # Define ticks do eixo x₂
  xticks!(fig, 0:(1/Nx1):1, rotation = 45)# Define ticks do eixo x₁
end
