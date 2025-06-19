struct Malha
  base # Informações sobre o tipo da base de funções interpoladoras
  ne::Int64 # Número de elementos totais
  neq::Int64 # Número de equações
  coords::Tuple # Tupla com as coordenadas (X₁, X₂, Xᵢ...) dos nós da malha
  dx # n-upla contendo o intervalo entre os nós, para cada eixo (se uniforme)
  EQ # Vetor com as reenumerações das funções globais ϕ
  LG # Matriz de conectividade local/global (LG). Relaciona a numeração local e global das funções ϕ
  EQoLG::Matrix{Int64} # Matriz composta entre a EQ e a LG
  a # Coordenada do início do intervalo uniforme
  b # Coordenada do final do intervalo uniforme
  n_dim # Número de dimensões da malha
  Nx # Número de subdivisões da malha para cada eixo
end

function montaLG_geral(Nx1::Int64, Nx2::Int64 = 0)
  n_dim = 1 + (Nx2 == 0 ? 0 : 1)
  ne = Nx1 * (Nx2 == 0 ? 1 : Nx2)

  LG = zeros(Int64, 2^n_dim, ne)

  # Monta o primeiro bloco
  LG1 = [i for i in 1:Nx1]
  for e in 1:Nx1
    for a in 1:(2^n_dim)
      @inbounds LG[a, e] = LG1[e] + ((a - 1) % 2) + (Nx1 + 1) * floor((a - 1) / 2)
    end
  end

  # Preenche o restante dos blocos baseado no primeiro
  for e in (Nx1 + 1):ne
    @inbounds LG[:, e] .= LG[:, e - Nx1] .+ (Nx1 + 1)
  end

  return LG
end

function montaLG(ne, base)
  LG = zeros(Int64, ne, base.nB)

  LG1 = [(base.nB - 1) * i - (base.nB - 2) for i in 1:ne]

  for e in 1:ne
    for a in 1:(base.nB)
      @inbounds LG[e, a] = LG1[e] + (a - 1)
    end
  end

  return LG'
end

function monta_LG_2D(Nx1::Int64, Nx2::Int64)::Matrix{Int64}
  # Define o número de funções φ no eixo x₁ e x₂
  nx1 = Nx1 + 1
  nx2 = Nx2 + 1

  # M[:,j] contém a numeração da primeira linha do "Bloco j"
  M = (1:(nx1 - 1)) .+ (0:nx1:((nx2 - 2) * nx1))'

  # LG[1,:] contém a numeração global da primeira função local de cada elemento 
  linha1 = reshape(M, 1, :)

  # Constrói a matriz LG
  LG = vcat(linha1, linha1 .+ 1, linha1 .+ nx1, linha1 .+ (nx1 + 1))

  return LG
end

function montaEQ_geral(Nx1::Int64, Nx2::Int64 = 0)
  n_dim = 1 + (Nx2 == 0 ? 0 : 1)
  ne = Nx1 * (Nx2 == 0 ? 1 : Nx2)

  neq = (Nx1 - 1) * (Nx2 == 0 ? 1 : Nx2 - 1)

  EQ = fill(neq + 1, (Nx1 + 1) * (Nx2 == 0 ? 1 : Nx2 + 1))

  contador = 1
  for y in 1:(Nx2 + 1)
    for x in 1:(Nx1 + 1)
      if !(x == 1 || x == Nx1 + 1 || (Nx2 > 0 && (y == 1 || y == Nx2 + 1)))
        idx = x + (y - 1) * (Nx1 + 1)
        EQ[idx] = contador
        contador += 1
      end
    end
  end

  return neq, EQ
end

function montaEQ(ne, neq, base)
  EQ = zeros(Int64, (base.nB - 1) * ne + 1)

  EQ[1] = neq + 1
  EQ[end] = neq + 1

  for i in 2:((base.nB - 1) * ne)
    @inbounds EQ[i] = i - 1
  end

  return EQ
end

function monta_EQ_2D(Nx1::Int64, Nx2::Int64)
  # Define o número de funções φ no eixo x₁ e x₂
  nx1 = Nx1 + 1
  nx2 = Nx2 + 1

  # Calcula o número de funções globais φ que compõem a base do espaço Vₘ
  m = (nx1 - 2) * (nx2 - 2)

  # Inicializa o vetor EQ preenchido com m+1
  EQ = fill(m + 1, nx1 * nx2)

  # Vetor contendo os índices das funções globais φ que compõem a base do espaço Vₘ
  L = reshape((0:(nx1 - 3)) .+ ((nx1 + 2):nx1:((nx2 - 2) * nx1 + 2))', :, 1)

  # Atribui os valores de 1 até m as funções globais φ que compõem a base do espaço Vₘ
  EQ[L] = 1:m

  return m, EQ
end

function monta_malha_1D_uniforme(ne, base, a, b)
  dx = (b[1] - a[1]) / ne
  X = collect(a:dx:b)

  neq, EQ = montaEQ_geral(ne)
  LG = montaLG_geral(ne)
  EQoLG = EQ[LG]

  coords = Tuple([X])

  n_dim = 1
  Nx = (; ne)
  return Malha(base, ne, neq, coords, dx, EQ, LG, EQoLG, a, b, n_dim, Nx)
end

function monta_malha_2D_uniforme(base, Nx1, Nx2, a::Tuple, b::Tuple)::Malha
  # Define o comprimento da base (h₁) e altura (h₂) de cada elemento retangular Ωᵉ
  h₁, h₂ = (b[1] - a[1]) / Nx1, (b[2] - a[2]) / Nx2
  h = (; h₁, h₂)

  ne = Nx1 * Nx2

  # Define a discretização em x₁ e x₂
  x₁ = collect(a[1]:h₁:b[1])
  x₂ = collect(a[2]:h₂:b[2])

  # Define as coordenadas de cada nó da malha
  X₁ = [x₁[i] for i in 1:(Nx1 + 1), j in 1:(Nx2 + 1)]
  X₂ = [x₂[j] for i in 1:(Nx1 + 1), j in 1:(Nx2 + 1)]

  neq, EQ = montaEQ_geral(Nx1, Nx2)
  LG = montaLG_geral(Nx1, Nx2)
  EQoLG = EQ[LG]

  coords = (X₁, X₂)
  n_dim = 2
  Nx = (; Nx1, Nx2)
  return Malha(base, ne, neq, coords, h, EQ, LG, EQoLG, a, b, n_dim, Nx)
end

function malha2D_adiciona_ruido(malha::Malha)::Malha
  (X₁, X₂) = malha.coords
  (; h₁, h₂) = malha.dx
  # Verifica se as matrizes têm a mesma dimensão
  @assert size(X₁)==size(X₂) "X₁ e X₂ devem ter as mesmas dimensões"

  noise_magnitude = 4
  # Verifica se a malha é suficientemente grande para ter nós internos
  if size(X₁, 1) > 2 && size(X₁, 2) > 2
    # Define os limites do ruído
    ruído_limite₁, ruído_limite₂ = h₁ / noise_magnitude, h₂ / noise_magnitude

    # Aplica ruído uniforme aos nós internos
    X₁[2:(end - 1),
    2:(end - 1)] .+= ruído_limite₁ *
                     (rand(Float64, size(X₁[2:(end - 1), 2:(end - 1)])) .- 0.5) * 2
    X₂[2:(end - 1),
    2:(end - 1)] .+= ruído_limite₂ *
                     (rand(Float64, size(X₂[2:(end - 1), 2:(end - 1)])) .- 0.5) * 2
  end

  return Malha(
    malha.base,
    malha.ne, malha.neq,
    (X₁, X₂), malha.dx,
    malha.EQ, malha.LG, malha.EQoLG,
    malha.a, malha.b,
    malha.n_dim,
    malha.Nx
  )
end