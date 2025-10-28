"""
    jacobian(n_dim::Int, Xᵉ::Vector, ∇Φξ::Tuple{Matrix{Float64}}, ξ::Int)

Calcula a matriz jacobiana, assim como seu determinante, 
no domínio do elemento determinado pelas coordenadas em `Xᵉ`.

# Parâmetros
- `n_dim::Int`: Número de dimensões do espaço.
- `Xᵉ::Vector`: Coordenadas dos vértices do elemento finito Ωᵉ. Cada entrada do vetor representa um eixo.
- `∇Φξ::Tuple{Matrix{Float64}}`: N-upla de matrizes `npg x n_funcs` com a avaliaçao das `∂ϕ` nos pontos de Gauss-Legendre.
- `ξ::Int`: Representa o número do ponto de Gauss.

# Retorno
- `M::Matrix{Float64}`: Matriz jacobiana.
- `detJ::Float64`: O determinante da matriz jacobiana

# Exemplo
```@example

```
"""
function jacobiano(n_dim, Xᵉ, ∇Φξ, ξ)
  # Monta a matriz jacobiana
  M = Matrix{Float64}(undef, n_dim, n_dim)
  for i in 1:n_dim
    for j in 1:n_dim
      M[i, j] = dot(Xᵉ[i], ∇Φξ[j][ξ, :])
    end
  end

  # Calcula o determinante Jacobiano
  detJ = det(M)

  return M, detJ
end

"""
    quadratura_gauss(npg::Int, n_dim::Int)

Define os `npg` pontos de Gauss-Legendre `P` e os pesos associados `W`. 
Cada um dos pontos em `P` são n-uplas dimensionadas por `n_dim`. Os pesos associados `W` funcionam analogamente.

# Parâmetros
- `npg::Int`: Número de pontos de integração de Gauss-Legendre.
- `n_dim::Int`: Número de dimensões do espaço. Tamanho dos pontos de Gauss-Legendre.

# Retorno
- `P::Vector{Tuple}`: Vetor de n-uplas, onde cada n-upla é um ponto de Gauss-Legendre dimensionado por `n_dim`.
- `W::Vector{Tuple}`: Vetor de n-uplas, onde cada n-upla é um conjunto ordenado de pesos de Gauss-Legendre associados a cada um dos pontos em `P`.

# Exemplo
```@example

```
"""
function quadratura_gauss(npg::Int, n_dim::Int)
  p, w = legendre(npg)

  # Declara as listas de pontos de Gauss e seus respectivos pesos
  P = Vector{Tuple}(undef, npg^n_dim)
  W = similar(P)

  # Os pontos de Gauss são definidos como combinações dos n_dim eixos
  # A ordem das combinações é do eixo de maior número para o menor
  # Exemplo para n_dim = 2: para trocar o ponto no eixo ξ₁ primeiro passa por todos os pontos no eixo ξ₂.
  for i in 1:(npg^n_dim)
    ponto = ()
    peso = ()
    for d in 1:n_dim
      idx_global = floor((Float64(npg)^(n_dim - d) + i - 1) * Float64(npg)^-(n_dim - d))
      idx = Int((idx_global - 1) % npg) + 1
      ponto = (ponto..., p[idx])
      peso = (peso..., w[idx])
    end
    P[i] = ponto
    W[i] = peso
  end

  return P, W
end

"""
    quadratura_ϕ(base, npg::Int, n_dim::Int)

Avalia cada uma das funções `ϕ` em cada um dos `npg` pontos de Gauss-Legendre e armazena o resultado em uma matriz `ϕP: npg x n_funcs`.

# Parâmetros
- `n_funcs::Int`: Número de funções base.
- `npg::Int`: Número de pontos de integração de Gauss-Legendre.
- `n_dim::Int`: Número de dimensões do espaço. Tamanho dos pontos de Gauss-Legendre.

# Retorno
- `ϕP::Matrix{Float64}`: Matriz `npg x n_funcs` com a avaliaçao das `ϕ` nos pontos de Gauss-Legendre.
- `P::Vector{Tuple}`: Vetor de n-uplas, onde cada n-upla é um ponto de Gauss-Legendre dimensionado por `n_dim`.
- `W::Vector{Tuple}`: Vetor de n-uplas, onde cada n-upla é um conjunto ordenado de pesos de Gauss-Legendre associados a cada um dos pontos em `P`.

# Exemplo
```@example

```
"""
function quadratura_ϕ(base, npg::Int, n_dim::Int) # Acho que esse base posso mudar para n_funcs
  P, W = quadratura_gauss(npg, n_dim)
  n_funcs = base.nB

  ϕP = zeros(Float64, npg^n_dim, n_funcs^n_dim)
  # Para todos os pontos de Gauss, avalia as ϕ locais
  for ξ in 1:(npg^n_dim)
    ϕP[ξ, :] .= ϕ_geral(P[ξ]...)[1]
  end

  return ϕP, P, W
end

"""
    quadratura_∇ϕ(base, npg::Int, n_dim::Int)

Avalia cada uma das funções `∂ϕ`, para cada uma das `n_dim` dimensões, em cada um dos `npg` pontos de Gauss-Legendre e 
armazena o resultado em uma n-upla `∇ϕP` de matrizes `∂ϕP: npg x n_funcs`.

# Parâmetros
- `n_funcs::Int`: Número de funções base.
- `npg::Int`: Número de pontos de integração de Gauss-Legendre.
- `n_dim::Int`: Número de dimensões do espaço. Tamanho dos pontos de Gauss-Legendre.

# Retorno
- `∇ϕP::Tuple{Matrix{Float64}}`: N-upla de matrizes `npg x n_funcs` com a avaliaçao das `∂ϕ` nos pontos de Gauss-Legendre.
- `P::Vector{Tuple}`: Vetor de n-uplas, onde cada n-upla é um ponto de Gauss-Legendre dimensionado por `n_dim`.
- `W::Vector{Tuple}`: Vetor de n-uplas, onde cada n-upla é um conjunto ordenado de pesos de Gauss-Legendre associados a cada um dos pontos em `P`.

# Exemplo
```@example

```
"""
function quadratura_∇ϕ(base, npg::Int, n_dim::Int) # Acho que esse base posso mudar para n_funcs
  P, W = quadratura_gauss(npg, n_dim)
  n_funcs = base.nB

  ∇ϕP = ()
  for d in 1:n_dim
    ∂ϕᵢP = zeros(npg^n_dim, n_funcs^n_dim)

    # Para todos os pontos de Gauss, avalia as ∂ϕᵢ locais (i = d)
    for ξ in 1:(npg^n_dim)
      ∂ϕᵢP[ξ, :] .= ∇ϕ_geral(P[ξ]...)[d]
    end
    ∇ϕP = (∇ϕP..., ∂ϕᵢP)
  end

  return ∇ϕP, P, W
end

"""
    mudanca_variavel_xξ(Xᵉ::Vector, Φξ::Vector, n_dim::Int)

Mapeia as coordenadas físicas x = (x₁, x₂, xᵢ...) a partir de ξ = (ξ₁, ξ₂, ξᵢ...)

# Parâmetros
- `Xᵉ::Vector`: Coordenadas dos vértices do elemento finito Ωᵉ. Cada entrada do vetor representa um eixo.
- `Φξ::Vector`: Vetor com o valor das funções locais Φ avaliadas no ponto de Gauss ξ = (ξ₁, ξ₂, ξᵢ...).
- `n_dim::Int`: Número de dimensões do espaço.

# Retorno
- `x::Vector`: A coordenada física referente ao ponto de Gauss ξ.

# Exemplo
```@example

```
"""
function mudanca_variavel_xξ(Xᵉ, Φξ, n_dim)
  x = ()

  # Mapeia as coordenadas físicas x = (x₁, x₂, xᵢ...) a partir de ξ = (ξ₁, ξ₂, ξᵢ...)
  for d in 1:n_dim
    x = (x..., dot(Xᵉ[d], Φξ))
  end
  return [x...]
end

"""
    elem_coords(malha::Malha, e::Int)

Encontra as coordenadas dos nós do elemento finito `e` e os retorna, junto a EQoLG referente.

# Parâmetros
- `malha::Malha`: A malha sob a qual está sendo solucionado o problema.
- `e::Int`: O identificador numérico do elemento finito Ωᵉ.

# Retorno
- `eqs_idx::`: EQoLG para o elemento `e` passado como parâmetro
- `Xᵉ::Vector`: Coordenadas dos vértices do elemento finito Ωᵉ. Cada entrada do vetor representa um eixo.

# Exemplo
```@example

```
"""
function elem_coords(malha::Malha, e::Int)
  (; LG, EQ, n_dim, coords) = malha

  # Índices dos vértices/funções globais do elemento finito Ωᵉ
  nosᵉ_idx = LG[:, e]

  # Coordenadas dos vértices do elemento finito Ωᵉ
  Xᵉ = []
  for d in 1:n_dim
    push!(Xᵉ, coords[d][nosᵉ_idx])
  end

  # Numeração de equação dos índices no vetor nosᵉ_idx
  eqs_idx = EQ[nosᵉ_idx]

  return eqs_idx, Xᵉ
end

"""
    montaKᵉ_geral!(Kᵉ, Xᵉ, P, W, Φξ, ∇Φξ, n_dim, pseudo_a)

Monta a matriz local `Kᵉ` referente à contribuição do elemento `e`.

# Parâmetros
- `Kᵉ::Matrix{Float64}`: Matriz local `Kᵉ` referente à contribuição do elemento `e`.
- `Xᵉ::Vector`: Coordenadas dos vértices do elemento finito Ωᵉ. Cada entrada do vetor representa um eixo.
- `P::Vector{Tuple}`: Vetor de n-uplas, onde cada n-upla é um ponto de Gauss-Legendre dimensionado por `n_dim`.
- `W::Vector{Tuple}`: Vetor de n-uplas, onde cada n-upla é um conjunto ordenado de pesos de Gauss-Legendre associados a cada um dos pontos em `P`.
- `Φξ::Matrix{Float64}`: Matriz `npg x n_funcs` com a avaliaçao das `ϕ` nos pontos de Gauss-Legendre.
- `∇Φξ::Tuple{Matrix{Float64}}`: N-upla de matrizes `npg x n_funcs` com a avaliaçao das `∂ϕ` nos pontos de Gauss-Legendre.
- `n_dim::Int`: Número de dimensões do espaço.
- `pseudo_a`: A referência ao operador bilinear a(u,v).

# Retorno
Altera `Kᵉ`.

# Exemplo
```@example

```
"""
function montaKᵉ_geral!(Kᵉ, Xᵉ, P, W, Φξ, ∇Φξ, n_dim, pseudo_a)
  # Zera as entradas da matriz local Kᵉ
  fill!(Kᵉ, 0.0)

  # Função que retorna o gradiente de Φ, dado um ponto de Gauss e a numeração da função local
  function ∇Φ(ξ, a)
    return [∇Φξ[dim][ξ, a] for dim in 1:n_dim]
  end

  npg = length(P)

  # Itera sobre os pontos de Gauss (combinações)
  for ξ in 1:npg
    # Vetor com o valor das funções locais Φ avaliadas no ponto de Gauss ξ = (ξ₁, ξ₂, ξᵢ...)
    ϕᵉ = Φξ[ξ, :]

    x = mudanca_variavel_xξ(Xᵉ, ϕᵉ, n_dim)

    # Matriz e determinante do Jacobiano
    M, detJ = jacobiano(n_dim, Xᵉ, ∇Φξ, ξ)
    @assert detJ>0 "O determinante jacobiano deve ser positivo"

    # Calcula a matriz H que aplica a mudança de variável de ∇ϕᵉ_a para ∇Φ_a
    M⁻¹ = inv(M)
    Hᵀ = M⁻¹ * detJ
    H = transpose(Hᵀ)

    # Calcula o peso total do ponto de Gauss ξ = (ξ₁, ξ₂, ξᵢ...)
    WW = prod(W[ξ])

    detJ⁻¹H = 1 / detJ * H

    # Calcula a contribuição de quadratura e acumula o valor na matriz Kᵉ
    @inbounds for a in 1:(2^n_dim)
      # Aplica a mudança de variável de ∇ϕᵉ_a para ∇Φ_a
      ∇ϕᵉ_a = detJ⁻¹H * ∇Φ(ξ, a)
      ϕᵉ_a = ϕᵉ[a]
      for b in 1:(2^n_dim)
        # Aplica a mudança de variável de ∇ϕᵉ_b para ∇Φ_b
        ∇ϕᵉ_b = detJ⁻¹H * ∇Φ(ξ, b)

        ϕᵉ_b = ϕᵉ[b]

        termos_equacao = TermosEquacao(
          ϕᵉ_b,
          ϕᵉ_a,
          ∇ϕᵉ_b,
          ∇ϕᵉ_a,
          x
        )
        soma = pseudo_a(termos_equacao)

        Kᵉ[a, b] += WW * soma * detJ
      end
    end
  end
end

"""
    montaK_geral(malha::Malha, pseudo_a)

Monta a matriz `K` global do sistema linear do problema.

# Parâmetros
- `malha::Malha`: A malha sob a qual está sendo solucionado o problema.
- `pseudo_a::`: A referência ao operador bilinear a(u,v).

# Retorno
- `K::Matrix{Float64}`: A matriz `K` global do sistema linear do problema.

# Exemplo
```@example

```
"""
function montaK_geral(malha::Malha, pseudo_a)
  (; ne, neq, n_dim, Nx, base) = malha

  npg = 2

  # Avalia as ϕ e ∇ϕ nos pontos de Gauss
  ϕξ, P, W = quadratura_ϕ(base, npg, n_dim)
  ∇ϕξ, P, W = quadratura_∇ϕ(base, npg, n_dim)

  # Define a banda da matriz K
  band = Nx[1]

  # Inicializa as matrizes locais e globais
  K = BandedMatrix(Zeros(neq, neq), (band, band))
  Kᵉ = zeros(Float64, 2^n_dim, 2^n_dim)

  # Itera sobre os elementos Ωᵉ
  for e in 1:ne
    eqs_idx, Xᵉ = elem_coords(malha::Malha, e::Int)

    # Calcula a matriz local Kᵉ
    montaKᵉ_geral!(Kᵉ, Xᵉ, P, W, ϕξ, ∇ϕξ, n_dim, pseudo_a)

    # Itera sobre as colunas (b) e linhas (a) da matriz local Kᵉ
    @inbounds for b in 1:(2^n_dim)
      j = eqs_idx[b]
      for a in 1:(2^n_dim)
        i = eqs_idx[a]
        if i <= neq && j <= neq
          K[i, j] += Kᵉ[a, b]
        end
      end
    end
  end

  return sparse(K)
end

"""
    montaFᵉ_geral!(Fᵉ, f, Xᵉ, P, W, ϕξ, ∇ϕξ, n_dim)

Monta o vetor local `Fᵉ` referente à contribuição do elemento `e`.

# Parâmetros
- `Fᵉ::Vector{Float64}`: Vetor local `Fᵉ` referente à contribuição do elemento `e`.
- `f::function`: Funçãp que represnenta o lado direito da equação do problema.
- `Xᵉ::Vector`: Coordenadas dos vértices do elemento finito Ωᵉ. Cada entrada do vetor representa um eixo.
- `P::Vector{Tuple}`: Vetor de n-uplas, onde cada n-upla é um ponto de Gauss-Legendre dimensionado por `n_dim`.
- `W::Vector{Tuple}`: Vetor de n-uplas, onde cada n-upla é um conjunto ordenado de pesos de Gauss-Legendre associados a cada um dos pontos em `P`.
- `ϕξ::Matrix{Float64}`: Matriz `npg x n_funcs` com a avaliaçao das `ϕ` nos pontos de Gauss-Legendre.
- `∇ϕξ::Tuple{Matrix{Float64}}`: N-upla de matrizes `npg x n_funcs` com a avaliaçao das `∂ϕ` nos pontos de Gauss-Legendre.
- `n_dim::Int`: Número de dimensões do espaço.

# Retorno
Altera `Fᵉ`.

# Exemplo
```@example

```
"""
function montaFᵉ_geral!(Fᵉ, f, Xᵉ, P, W, ϕξ, ∇ϕξ, n_dim)
  # Zera as entradas do vetor local Fᵉ
  fill!(Fᵉ, 0.0)

  npg = length(P)

  # Itera sobre os pontos de Gauss (combinações)
  for ξ in 1:npg
    # Vetor com o valor das funções locais Φ avaliadas no ponto de Gauss ξ = (ξ₁, ξ₂, ξᵢ...)
    ϕᵉ = ϕξ[ξ, :]

    x = mudanca_variavel_xξ(Xᵉ, ϕᵉ, n_dim)

    # Matriz e determinante do Jacobiano
    M, detJ = jacobiano(n_dim, Xᵉ, ∇ϕξ, ξ)
    @assert detJ>0 "O determinante jacobiano deve ser positivo"

    # Calcula o peso total do ponto de Gauss ξ = (ξ₁, ξ₂, ξᵢ...)
    WW = prod(W[ξ])

    # Calcula a contribuição de quadratura e acumula o valor no vetor Fᵉ
    for a in 1:(2^n_dim)
      @inbounds Fᵉ[a] += WW * f(x...) * ϕᵉ[a] * detJ
    end
  end
end

"""
    montaF_geral(f::Function, malha::Malha)

Monta o vetor `F` global do sistema linear do problema.

# Parâmetros
- `f::function`: Função que represnenta o lado direito da equação do problema.
- `malha::Malha`: A malha sob a qual está sendo solucionado o problema.

# Retorno
- `F::Vector{Float64}`: O vetor `F` global do sistema linear do problema.

# Exemplo
```@example

```
"""
function montaF_geral(f::Function, malha::Malha)
  (; ne, neq, n_dim, base) = malha

  npg = 5

  # Avalia as ϕ e ∇ϕ nos pontos de Gauss
  ϕξ, P, W = quadratura_ϕ(base, npg, n_dim)
  ∇ϕξ, P, W = quadratura_∇ϕ(base, npg, n_dim)

  # Inicializa os vetores locais e globais
  F = zeros(neq + 1)
  Fᵉ = zeros(2^n_dim)

  # Itera sobre os elementos Ωᵉ
  for e in 1:ne
    eqs_idx, Xᵉ = elem_coords(malha::Malha, e::Int)

    # Calcula o vetor local Fᵉ
    montaFᵉ_geral!(Fᵉ, f, Xᵉ, P, W, ϕξ, ∇ϕξ, n_dim)

    # Adiciona a contribuição do elemento finito ao vetor global F
    for a in 1:(2^n_dim)
      i = eqs_idx[a]
      F[i] += Fᵉ[a]
    end
  end

  return F[1:neq]
end

"""
    solve_sys_poisson(run_values::RunValues, malha::Malha)

Monta e soluciona o sistema linear KC = F, considerando o problema com equação de Poisson.

# Parâmetros
- `run_values::RunValues`: Estrutura que contém parâmetros do problema.
- `malha::Malha`: A malha sob a qual está sendo solucionado o problema.

# Retorno
- `C::Vector{Float64}`: O vetor contendo a solução do sistema linear.

# Exemplo
```@example

```
"""
function solve_sys_poisson(run_values::RunValues, malha::Malha)
  (; α, β, f) = run_values

  function pseudo_a(termos_equacao::TermosEquacao)
    (; ∇u, ∇v, u, v) = termos_equacao

    return β * dot(u, v) + α * dot(∇u, ∇v)
  end

  C = solve_sys(f, malha, pseudo_a)

  return C
end

"""
    solve_sys(f::function, malha::Malha, pseudo_a::function)

Monta e soluciona o sistema linear KC = F, 
considerando a referência ao operador bilinear a(u,v) obtido na formulação fraca.

# Parâmetros
- `run_values::RunValues`: Estrutura que contém parâmetros do problema.
- `malha::Malha`: A malha sob a qual está sendo solucionado o problema.

# Retorno
- `C::Vector{Float64}`: O vetor contendo a solução do sistema linear.

# Exemplo
```@example

```
"""
function solve_sys(f, malha, pseudo_a)
  K = montaK_geral(malha, pseudo_a)

  F = montaF_geral(f, malha)

  C = zeros(Float64, malha.neq)

  C .= K \ F

  return C
end

"""
    monta_u_aproximada(c::Vector{Float64}, malha::Malha)

Dada a solução obtida através do sistema linear e a malha sob a qual foi solucionada,
retorna uma função que aproxima a solução analítica do problema.

# Parâmetros
- `c::Vector{Float64}`: O vetor contendo a solução do sistema linear.
- `malha::Malha`: A malha sob a qual está sendo solucionado o problema.

# Retorno
- `uh::function`: Uma função que aproxima a solução analítica do problema.

# Exemplo
```@example

```
"""
function monta_u_aproximada(c, malha::Malha)
  (; ne, n_dim) = malha
  d = [c..., 0]

  ξₓ(x, x₀, x_f) = 2 * (x - x₀) / (x_f - x₀) - 1
  function uh(x...)
    soma = 0
    x₀ = Array{Float64}(undef, 0)
    x_f = Array{Float64}(undef, 0)

    for e in 1:ne
      j, Xᵉ = elem_coords(malha::Malha, e::Int)
      for dim in 1:n_dim
        X = Xᵉ[dim]
        push!(x₀, X[1])
        push!(x_f, X[end])
      end

      if all((x .- x₀) .>= 0) && all((x_f .- x) .> 0)
        ξ = ξₓ.(x, x₀, x_f)

        soma += dot(d[j], ϕ_geral(ξ...)[1])
      end

      empty!(x₀)
      empty!(x_f)
    end

    return soma
  end

  return uh
end
