
"""
    struct Mesh

Valores que definem uma malha.

# Campos
- `base::Base`: Informações sobre o tipo da base de funções interpoladoras
- `ne::Int64`: Número de elementos totais
- `neq::Int64`: Número de equações
- `coords::Tuple`: N-upla com as coordenadas (X₁, X₂, Xᵢ...) dos nós da malha
- `dx::`: n-upla contendo o intervalo entre os nós, para cada eixo (se uniforme)
- `EQ::Vector{Int64}`: Vetor com as reenumerações das funções globais ϕ
- `LG::Matrix{Int64}`: Matriz de conectividade local/global (LG). Relaciona a numeração local e global das funções ϕ
- `EQoLG::Matrix{Int64}`: Matriz composta entre a EQ e a LG
- `a::Tuple`: Coordenada do início do intervalo
- `b::Tuple`: Coordenada do final do intervalo
- `n_dim::Int64`: Número de dimensões da malha
- `Nx::Tuple`: Número de subdivisões da malha para cada eixo
"""
struct Malha
  base
  ne::Int64
  neq::Int64
  coords::Tuple
  dx
  EQ::Vector{Int64}
  LG::Matrix{Int64}
  EQoLG::Matrix{Int64}
  a
  b
  n_dim::Int64
  Nx
end

"""
    montaLG_geral(Nx1::Int64, Nx2::Int64 = 0)

Função que monta a matriz de conectividade local/global (LG). Relaciona a numeração local e global das funções ϕ.

# Parâmetros
- `Nx1::Int64`: Número de subdivisões da malha no primeiro eixo.
- `Nx2::Int64 = 0`: Número de subdivisões da malha no segundo eixo.

# Retorno
- `LG::Matrix{Int64}`: Matriz de conectividade local/global (LG). Relaciona a numeração local e global das funções ϕ.

# Exemplo
```@example

```
"""
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
  for e in (Nx1+1):ne
    @inbounds LG[:, e] .= LG[:, e-Nx1] .+ (Nx1 + 1)
  end

  return LG
end

"""
    monta_LG_2D(Nx1::Int64, Nx2::Int64)::Matrix{Int64}

Função que monta a matriz de conectividade local/global (LG). Relaciona a numeração local e global das funções ϕ.

# Parâmetros
- `Nx1::Int64`: Número de subdivisões da malha no primeiro eixo.
- `Nx2::Int64`: Número de subdivisões da malha no segundo eixo.

# Retorno
- `LG::Matrix{Int64}`: Matriz de conectividade local/global (LG). Relaciona a numeração local e global das funções ϕ.

# Exemplo
```@example

```
"""
function monta_LG_2D(Nx1::Int64, Nx2::Int64)::Matrix{Int64}
  # Define o número de funções φ no eixo x₁ e x₂
  nx1 = Nx1 + 1
  nx2 = Nx2 + 1

  # M[:,j] contém a numeração da primeira linha do "Bloco j"
  M = (1:(nx1-1)) .+ (0:nx1:((nx2-2)*nx1))'

  # LG[1,:] contém a numeração global da primeira função local de cada elemento 
  linha1 = reshape(M, 1, :)

  # Constrói a matriz LG
  LG = vcat(linha1, linha1 .+ 1, linha1 .+ nx1, linha1 .+ (nx1 + 1))

  return LG
end

"""
    montaEQ_geral(Nx1::Int64, Nx2::Int64 = 0)

Função que monta o vetor com as reenumerações das funções globais ϕ.

# Parâmetros
- `Nx1::Int64`: Número de subdivisões da malha no primeiro eixo.
- `Nx2::Int64 = 0`: Número de subdivisões da malha no segundo eixo.

# Retorno
- `EQ::Vector{Int64}`: Vetor com as reenumerações das funções globais ϕ.

# Exemplo
```@example

```
"""
function montaEQ_geral(Nx1::Int64, Nx2::Int64 = 0)
  n_dim = 1 + (Nx2 == 0 ? 0 : 1)
  ne = Nx1 * (Nx2 == 0 ? 1 : Nx2)

  neq = (Nx1 - 1) * (Nx2 == 0 ? 1 : Nx2 - 1)

  EQ = fill(neq + 1, (Nx1 + 1) * (Nx2 == 0 ? 1 : Nx2 + 1))

  contador = 1
  for y in 1:(Nx2+1)
    for x in 1:(Nx1+1)
      if !(x == 1 || x == Nx1 + 1 || (Nx2 > 0 && (y == 1 || y == Nx2 + 1)))
        idx = x + (y - 1) * (Nx1 + 1)
        EQ[idx] = contador
        contador += 1
      end
    end
  end

  return neq, EQ
end

"""
    monta_malha_1D_uniforme(baseType, Nx1, a, b)::Malha

Função que constrói uma malha 1D uniforme.

# Parâmetros
- `baseType::`: O tipo da base das funções interpoladoras.
- `Nx1::Int64`: Número de subdivisões da malha no primeiro eixo.
- `a::Tuple`: Coordenada do início do intervalo uniforme.
- `b::Tuple`: Coordenada do final do intervalo uniforme.

# Retorno
- `malha::Malha`: Uma malha 1D uniforme dentro do intervalo `[a,b]`.

# Exemplo
```@example

```
"""
function monta_malha_1D_uniforme(baseType, ne, a, b)
  dx = (b[1] - a[1]) / ne
  X = collect(a:dx:b)

  neq, EQ = montaEQ_geral(ne)
  LG = montaLG_geral(ne)
  EQoLG = EQ[LG]

  coords = Tuple([X])

  n_dim = 1
  Nx = (; ne)

  base = monta_base(baseType, ne)
  return Malha(base, ne, neq, coords, dx, EQ, LG, EQoLG, a, b, n_dim, Nx)
end

"""
    monta_malha_2D_uniforme(baseType, Nx1, Nx2, a::Tuple, b::Tuple)::Malha

Função que constrói uma malha 2D uniforme.

# Parâmetros
- `baseType::`: O tipo da base das funções interpoladoras.
- `Nx1::Int64`: Número de subdivisões da malha no primeiro eixo.
- `Nx2::Int64`: Número de subdivisões da malha no segundo eixo. 
- `a::Tuple`: Coordenada do início do intervalo
- `b::Tuple`: Coordenada do final do intervalo 

# Retorno
- `malha::Malha`: Uma malha 2D uniforme dentro do intervalo `[a,b]`.

# Exemplo
```@example

```
"""
function monta_malha_2D_uniforme(baseType, Nx1, Nx2, a::Tuple, b::Tuple)::Malha
  # Define o comprimento da base (h₁) e altura (h₂) de cada elemento retangular Ωᵉ
  h₁, h₂ = (b[1] - a[1]) / Nx1, (b[2] - a[2]) / Nx2
  h = (; h₁, h₂)

  ne = Nx1 * Nx2

  # Define a discretização em x₁ e x₂
  x₁ = collect(a[1]:h₁:b[1])
  x₂ = collect(a[2]:h₂:b[2])

  # Define as coordenadas de cada nó da malha
  X₁ = [x₁[i] for i in 1:(Nx1+1), j in 1:(Nx2+1)]
  X₂ = [x₂[j] for i in 1:(Nx1+1), j in 1:(Nx2+1)]

  neq, EQ = montaEQ_geral(Nx1, Nx2)
  LG = montaLG_geral(Nx1, Nx2)
  EQoLG = EQ[LG]

  coords = (X₁, X₂)
  n_dim = 2
  Nx = (; Nx1, Nx2)

  base = monta_base(baseType, ne)
  return Malha(base, ne, neq, coords, h, EQ, LG, EQoLG, a, b, n_dim, Nx)
end

"""
    monta_malha_2D_uniforme(baseType, Nx1, Nx2, a::Tuple, b::Tuple)::Malha

Função que constrói uma malha 2D uniforme.

# Parâmetros
- `baseType::`: O tipo da base das funções interpoladoras.
- `Nx1::Int64`: Número de subdivisões da malha no primeiro eixo.
- `Nx2::Int64`: Número de subdivisões da malha no segundo eixo. 
- `a::Tuple`: Coordenada do início do intervalo
- `b::Tuple`: Coordenada do final do intervalo 

# Retorno
- `malha::Malha`: Uma malha 2D uniforme dentro do intervalo `[a,b]`.

# Exemplo
```@example

```
"""
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
    X₁[2:(end-1),
    2:(end-1)] .+= ruído_limite₁ *
                   (rand(Float64, size(X₁[2:(end-1), 2:(end-1)])) .- 0.5) * 2
    X₂[2:(end-1),
    2:(end-1)] .+= ruído_limite₂ *
                   (rand(Float64, size(X₂[2:(end-1), 2:(end-1)])) .- 0.5) * 2
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

"""
    monta_malha_2D_foco(
    baseType, Nx1, Nx2, a::Tuple, b::Tuple, ponto_foco::Tuple, precisao::Int64)::Malha

Função que constrói uma malha 2D que cria mais elementos próximos a um ponto especificado.

# Parâmetros
- `baseType::`: O tipo da base das funções interpoladoras.
- `Nx1::Int64`: Número de subdivisões da malha no primeiro eixo.
- `Nx2::Int64`: Número de subdivisões da malha no segundo eixo. 
- `a::Tuple`: Coordenada do início do intervalo
- `b::Tuple`: Coordenada do final do intervalo 
- `ponto_foco::Tuple`: O ponto na malha onde elementos serão adicionados em volta.
- `precisao::Int64 = 1`: A quantidade de iterações de incremento de intervalos em torno do ponto_foco

# Retorno
- `malha::Malha`: Uma malha 2D com concentração de elementos em `ponto_foco`.

# Exemplo
```@example

```
"""
function monta_malha_2D_foco(
    baseType, Nx1, Nx2, a::Tuple, b::Tuple, ponto_foco::Tuple, precisao::Int64 = 1)::Malha

  # Define o comprimento da base (h₁) e altura (h₂) de cada elemento retangular Ωᵉ
  h₁, h₂ = (b[1] - a[1]) / Nx1, (b[2] - a[2]) / Nx2
  h = (; h₁, h₂)

  # Define a discretização em x₁ e x₂
  x₁ = collect(a[1]:h₁:b[1])
  x₂ = collect(a[2]:h₂:b[2])

  (p₁, p₂) = ponto_foco

  # Encontra nó central para ambos eixos
  idx₁ = findfirst(x -> x == p₁, x₁)
  idx₂ = findfirst(x -> x == p₂, x₂)

  # Se não existir nó central, em algum dois eixos, cria-o
  if idx₁ == nothing
    idx₁ = findfirst(x -> x > p₁, x₁)
    insert!(x₁, idx₁, p₁)
    Nx1 += 1
  end
  if idx₂ == nothing
    idx₂ = findfirst(x -> x > p₂, x₂)
    insert!(x₂, idx₂, p₂)
    Nx2 += 1
  end

  # Para cada eixo, adiciona um nó antes e depois do nó central.
  # Esses novos nós se encontram na média entre o nó central e seus vizinhos.
  for i in 1:precisao
    # Encontra as médias entre nó central e vizinhos
    prev₁ = (x₁[idx₁-1] + x₁[idx₁]) / 2
    next₁ = (x₁[idx₁+1] + x₁[idx₁]) / 2

    # Insere nós entre vizinhos e central
    insert!(x₁, idx₁, prev₁)
    idx₁ += 1 # Ajusta valor do indice do nó central
    insert!(x₁, idx₁ + 1, next₁)
    Nx1 += 2 # Ajusta valor da quantidade de intervalos no eixo x₁

    # Análogo ao trecho acima
    prev₂ = (x₂[idx₂-1] + x₂[idx₂]) / 2
    next₂ = (x₂[idx₂+1] + x₂[idx₂]) / 2
    insert!(x₂, idx₂, prev₂)
    idx₂ += 1
    insert!(x₂, idx₂ + 1, next₂)
    Nx2 += 2
  end

  # Define as coordenadas de cada nó da malha
  X₁ = [x₁[i] for i in 1:(Nx1+1), j in 1:(Nx2+1)]
  X₂ = [x₂[j] for i in 1:(Nx1+1), j in 1:(Nx2+1)]

  neq, EQ = montaEQ_geral(Nx1, Nx2)
  LG = montaLG_geral(Nx1, Nx2)
  EQoLG = EQ[LG]

  coords = (X₁, X₂)
  n_dim = 2
  Nx = (; Nx1, Nx2)

  ne = Nx1 * Nx2
  base = monta_base(baseType, ne)
  return Malha(base, ne, neq, coords, h, EQ, LG, EQoLG, a, b, n_dim, Nx)
end
