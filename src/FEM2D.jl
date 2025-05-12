function exemplo1()
	alpha = 1.0
	beta = 1.0
	f(x₁,x₂) = (2*alpha*π^2+beta) * sin(π*x₁) * sin(π*x₂)
	u(x₁,x₂) = sin(π*x₁) * sin(π*x₂)

	return alpha, beta, f, u
end

############### TESTES ###############
function plot_solução_aproximada(c̄::Vector{Float64}, Nx1::Int64, Nx2::Int64, EQoLG::Matrix{Int64})
  # Comprimentos da base (h₁) e altura (h₂) de cada elemento retangular Ωᵉ
  h₁, h₂ = 1 / Nx1, 1 / Nx2

  # Define uma discretização do intervalo de referência [-1, 1] nos eixos ξ₁ e ξ₂
  ξ₁ = collect(-1:0.1:1) 
  ξ₂ = collect(-1:0.1:1) 

  # Inicializa um objeto de gráfico que acumulará as superfícies
  plt = Plots.plot(seriestype = :surface, title="Solução Aproximada",
  color=:viridis,
  xlabel="x₁", ylabel="x₂", colorbar=false, zlims=(0,1),
  xticks=[0, 0.5, 1], yticks=[0, 0.5, 1], zticks=[0, 1])

  # Loop nos elementos da malha (percorrendo cada subdivisão ao longo de x₂ e x₁)
  for j = 1:Nx2
    # Define a segunda coordenada do ponto inferior esquerdo do retângulo `Ωᵉ`.
    p₂ = (j - 1) * h₂
    # Calcula os valores correspondentes no eixo x₂ para os pontos ξ₂ em [-1,1]
    x₂ = x₂_de_ξ.(ξ₂, h₂, p₂)

    for i = 1:Nx1
      # Primeira coordenada do ponto inferior esquerdo do retângulo `Ωᵉ`.
      p₁ = (i - 1) * h₁
      # Valores correspondentes no eixo x₁ para os pontos ξ₁ em [-1,1]
      x₁ = x₁_de_ξ.(ξ₁, h₁, p₁)

      # Determina a numeração do elemento finito atual (`e`) na malha
      e = (j - 1) * Nx1 + i

      # Obtém os coeficientes `c` da solução aproximada no elemento `e`
      ce = c̄[EQoLG[:, e]]

      # Solução aproximada na malha x₁ × x₂. Tamanho length(x₂) x length(x₁)
      mat = [dot(ce,ϕ_2D(ξ₁[a],ξ₂[b])) for b in 1:length(x₂), a in 1:length(x₁)]

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
      Plots.plot!(plt, x₁, x₂, mat, seriestype = :surface, alpha=0.5)
    end
  end

  return plt
end

function malha2D_adiciona_ruido!(X₁::Matrix{Float64}, X₂::Matrix{Float64}, h₁::Float64, h₂::Float64)
  # Verifica se as matrizes têm a mesma dimensão
  @assert size(X₁) == size(X₂) "X₁ e X₂ devem ter as mesmas dimensões"
  
  # Verifica se a malha é suficientemente grande para ter nós internos
  if size(X₁, 1) > 2 && size(X₁, 2) > 2
    # Define os limites do ruído
    ruído_limite₁, ruído_limite₂ = h₁ / 4, h₂ / 4
    
    # Aplica ruído uniforme aos nós internos
    X₁[2:end-1, 2:end-1] .+= ruído_limite₁ * 
              (rand(Float64, size(X₁[2:end-1,2:end-1])) .- 0.5) * 2
    X₂[2:end-1, 2:end-1] .+= ruído_limite₂ * 
                        (rand(Float64, size(X₂[2:end-1,2:end-1])) .- 0.5) * 2
  end
end

function plot_malha2D(X₁::AbstractArray{Float64}, X₂::AbstractArray{Float64}, LG::Matrix{Int64})
  # Cria uma nova figura
  fig = Plots.plot( legend=false, aspect_ratio=:equal,
                    xticks=0:0.25:1,  # Defina os ticks do eixo x
                    yticks=0:0.25:1)  # Defina os ticks do eixo y

  # Adiciona os nós da malha
  Plots.scatter!(X₁, X₂, markersize=4, color=:blue)

  # Adiciona as linhas de contorno dos elementos
  for e in 1:size(LG, 2)
    # Obtém os índices dos nós do elemento `e`
    i₁, i₂, i₃, i₄ = LG[:, e]

    # Adiciona as linhas de contorno do elemento `e`
    Plots.plot!([X₁[i₁], X₁[i₂], X₁[i₄], X₁[i₃], X₁[i₁]], 
                [X₂[i₁], X₂[i₂], X₂[i₄], X₂[i₃], X₂[i₁]], color=:black)
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
  Nx1 = 4; Nx2 = 3;

	fig = exemplo_malha2D(Nx1,Nx2; ruido = false)  # Aqui a função é chamada inicialmente
	xlabel!(fig,"x_1")
	ylabel!(fig,"x_2")
	yticks!(fig, 0:1/Nx2:1)             # Define ticks do eixo x₂
	xticks!(fig, 0:1/Nx1:1, rotation=45)# Define ticks do eixo x₁
end