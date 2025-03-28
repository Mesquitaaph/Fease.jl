using LinearAlgebra, Plots, GaussQuadrature, DataFrames, SparseArrays

function exemplo1()
  α = 1.0
  β = 1.0
  f(x) = -2.0 * α + β * x *(x-1.0)
  u(x) = x * (x - 1.0)
  return α, β, f, u
end
function Exemplo2()
  α = 1.0 
  β = 1.0
  f(x) = (α*π^2 + β) *sin(π*x)
  u(x) = sin(π*x)
  return α, β, f, u
end
function Exemplo3()
  α = 1.0 
  β = 0.0
  f(x) = 8.0
  u(x) = -4.0 * x * (x - 1.0)
  return α, β, f, u
end
function Exemplo4()
  α = 1.0 
  β = 1.0
  f(x) = x
  u(x) = x + (exp(-x) - exp(x)) / (exp(1) - exp(-1))
  return α, β, f, u
end

"Retorna os valores de `α`, `β`, `f` e `u` de acordo com o exemplo escolhido."
function dados_de_entrada(exemplo::Int)
  if exemplo == 1
      return exemplo1()
  elseif exemplo == 2
      return Exemplo2()
  elseif exemplo == 3
      return Exemplo3()
  elseif exemplo == 4
      return Exemplo4()
  else
      error("Exemplo inválido. Escolha 1, 2, 3 ou 4.")
  end
end

"Retorna o valor de `ϕₐ(ξ)` para um dado `ξ` e `a`."
function ϕ(ξ::Float64, a::Int) :: Float64
    if a == 1
        return (1.0 - ξ) / 2.0
    elseif a == 2
        return (1.0 + ξ) / 2.0
    else
        error("a deve ser 1 ou 2.")
    end
end

"Retorna a derivada de `ϕₐ(ξ)` em relação a `ξ` para um dado `ξ` e `a`."
function dϕ(ξ, a)
    if a == 1
        return -1.0 / 2.0
    elseif a == 2
        return 1.0 / 2.0
    else
        error("a deve ser 1 ou 2.")
    end
end

"Retorna a matriz local `Kᵉ` de um elemento arbitrário `e` de tamanho `h`."
function monta_Kᵉ(α::Float64, β::Float64, h::Float64) :: Matrix{Float64}
    # Constantes associadas aos termos da primeira e da segunda integral
    cst1 = 2.0 * α / h
    cst2 = β * h / 2.0

    # Pontos e pesos de quadratura de Gauss-Legendre
    n_pg = 2
    P, W = legendre(n_pg)

    # Inicializa a matriz local Kᵉ
    Kᵉ = zeros(Float64, 2, 2)

    # Loop sobre as entradas da matriz local Kᵉ
    for a = 1:2
        for b = 1:2
            # Loop sobre os pontos de quadratura
            for j = 1:n_pg
                # Calcula a contribuição da quadratura na entrada Kᵉ[a,b]
                Kᵉ[a,b] += W[j] * ( cst1 * dϕ(P[j],a) * dϕ(P[j],b) + cst2 * ϕ(P[j],a) * ϕ(P[j],b) )
            end
        end
    end

    return Kᵉ
end

"""
Mapeia um ponto `ξ` do intervalo padrão `[-1, 1]` para o intervalo físico `[x1e, x2e]`.

# Parâmetros
- `ξ::Float64`: Ponto no intervalo padrão `[-1, 1]`.
- `h::Float64`: Comprimento do e-ésimo elemento, calculado como `h = x2e - x1e`.
- `x1e::Float64`: Coordenada do ponto inicial `x1e` do e-ésimo elemento.

# Retorno
- Coordenada correspondente a `ξ` no intervalo `[x1e, x2e]`.
"""
function x_de_ξ(ξ::Float64, h::Float64, x1e::Float64) :: Float64
    return (ξ + 1.0) * h / 2.0 + x1e
end

"""
Monta o vetor força local `Fᵉ` para o `e`-ésimo elemento `[x1e, x2e]` de tamanho `h`.

# Parâmetros
- `f::Function`: Função que define a força aplicada no intervalo.
- `h::Float64`: Comprimento do e-ésimo elemento, calculado como `h = x2e - x1e`.
- `x1e::Float64`: Coordenada do ponto inicial `x1e` do e-ésimo elemento.
- `P::Vector{Float64}`: Vetor dos pontos de quadratura no intervalo padrão `[-1, 1]`.
- `W::Vector{Float64}`: Vetor dos pesos de quadratura associados a `P`.

# Retorno
- `Fᵉ::Vector{Float64}`: Vetor força local de tamanho 2, correspondente ao elemento `[x1e, x2e]`.
"""
function monta_Fᵉ(f::Function, h::Float64, x1e::Float64, P::Vector{Float64}, W::Vector{Float64}) :: Vector{Float64}
    # Inicializa o vetor local Fᵉ
    Fᵉ = zeros(Float64, 2)

    # Loop sobre as entradas do vetor local Fᵉ
    for a = 1:2
        # Loop sobre os pontos de quadratura
        for j=1:length(P)
            # Mapeia o ponto de quadratura `P[j]` de [-1,1] para [x1e,x2e]
            x = x_de_ξ(P[j], h, x1e)
            # Calcula e acumula a contribuição da quadratura em Fᵉ[a]
            Fᵉ[a] += W[j] * f(x) * ϕ(P[j],a) * h/2
        end
    end

    return Fᵉ
end

"""
Retorna uma matriz `2 x ne` com a numeração global de cada função local.

# Parâmetros
- `ne::Int`: Número total de elementos (ou intervalos) em que o domínio `[0,1]` é dividido.
"""
function monta_LG(ne::Int) :: Matrix{Int}
    return transpose(hcat(1:ne, 2:ne+1))
end

"""
Retorna o valor de `m` o vetor `EQ`.

# Parâmetros
- `ne::Int`: Número total de elementos (ou intervalos) em que o domínio `[0,1]` é dividido.

# Retorno
- Uma tupla contendo:
  - `m::Int`: A dimensão do espaço aproximado `Vₘ`. Neste problema temos que `m = ne - 1`.
  - `EQ::Vector{Int}`: Um vetor de inteiros que fornece a reenumeração das funções globais `φ`. Neste problema temos que EQ = [m + 1; 1:m; m + 1].
"""
function monta_EQ(ne::Int) :: Tuple{Int,Vector{Int}}
    m = ne - 1         
    EQ = vcat(m+1,1:m,m+1)
    return m, EQ
end

"""
Retorna a matriz tridiagonal simétrica `K` de tamanho `m x m` montada a partir da matriz local `Kᵉ`.

# Parâmetros
- `α::Float64`: Parâmetro α fornecido como dado de entrada da equação diferencial em estudo.
- `β::Float64`: Parâmetro β fornecido como dado de entrada da equação diferencial em estudo.
- `ne::Int`: Número total de elementos (ou intervalos) em que o domínio `[0,1]` é dividido.
- `m::Int`: Dimensão do espaço aproximado `Vₘ`.
- `EQoLG::Matrix{Int}`: ≡ EQ[LG].
"""
function monta_K(α::Float64, β::Float64, ne::Int, m::Int, EQoLG::Matrix{Int}) :: SparseMatrixCSC{Float64, Int64}
    h = 1.0 / ne           # Comprimento de cada elemento finito
    K = spzeros(m+1, m+1)  # Inicializa a matriz esparsa K com tamanho (m+1) x (m+1)
    Kᵉ = monta_Kᵉ(α, β, h) # Calcula a matriz local Kᵉ

    # Loop sobre os elementos
    for e = 1:ne
        for a = 1:2
            i = EQoLG[a,e]
            for b = 1:2
                j = EQoLG[b,e]
                K[i,j] += Kᵉ[a,b]
            end
        end
    end
   
    # Remove a última linha e coluna da matriz K
    return K[1:m, 1:m]
end

"""
Retorna o vetor força `F` de tamanho `m`, utilizando os vetores locais `Fᵉ`.

# Parâmetros
- `f::Function`: Função que define a força aplicada no intervalo.
- `ne::Int`: Número total de elementos (ou intervalos) em que o domínio `[0,1]` é dividido.
- `m::Int`: Dimensão do espaço aproximado `Vₘ`.
- `EQoLG::Matrix{Int}`: ≡ EQ[LG].
"""
function monta_F(f::Function, ne::Int, m::Int, EQoLG::Matrix{Int}) :: Vector{Float64}
    h = 1.0 / ne       # Comprimento de cada elemento finito
    F = zeros(m+1)     # Inicializa o vetor F
    P, W = legendre(5) # Pontos e pesos de quadratura de Gauss-Legendre
    
    # Loop sobre os elementos
    for e = 1:ne
        Fᵉ = monta_Fᵉ(f, h, (e-1)*h, P, W) # Calcula o vetor local Fᵉ
        for a = 1:2
            i = EQoLG[a,e]
            F[i] += Fᵉ[a]
        end
    end

    # Remove a última entrada do vetor global F
    return F[1:m]
end

@doc raw"""
Calcula o erro na norma do espaço ``L^2`` entre a solução exata ``u(x)`` e a solução aproximada ``u_h(x)``, 
definida por ``\displaystyle u_h(x) = \sum_{j=1}^m c_j \varphi_j(x)``, em um domínio particionado em ``n_e`` elementos finitos. 

# Parâmetros:
- `u::Function`: Função que representa a solução exata ``u(x)``.
- `c̄::Vector{Float64}`: Vetor com os coeficientes da solução aproximada acrescido de um zero, i.e., `c̄ = [c; 0]`.
- `ne::Int`: Número total de elementos (ou intervalos) em que o domínio `[0,1]` é dividido.
- `EQoLG::Matrix{Int}`: ≡ EQ[LG].

# Retorno:
- `Float64`: O erro calculado na norma ``L^2``.
"""
function erro_norma_L2(u::Function, c̄::Vector{Float64}, ne::Int, EQoLG::Matrix{Int}) :: Float64
    h = 1.0 / ne       # Comprimento de cada elemento finito
    erro = 0.0         # Inicializa o erro
    P, W = legendre(5) # Pontos e pesos da quadratura de Gauss-Legendre (5 pontos)

    # Calcula as funções de base ϕ₁ e ϕ₂ nos pontos de quadratura P
    ϕ1P = ϕ.(P, 1)
    ϕ2P = ϕ.(P, 2)

    # Loop sobre todos os elementos finitos
    for e = 1:ne
        # Coeficientes da solução aproximada para o elemento `e`
        c1e = c̄[EQoLG[1,e]]
        c2e = c̄[EQoLG[2,e]]

        # Loop sobre os pontos de quadratura
        for j = 1:length(P)
            x = x_de_ξ(P[j], h, (e-1)*h) # Mapeia o ponto de quadratura `P[j]` de [-1,1] para [x1e,x2e]
            erro += W[j] * (u(x) - c1e*ϕ1P[j] - c2e*ϕ2P[j])^2
        end
    end

    return sqrt(erro*h/2)
end

# Carrega o exemplo escolhido
exemplo = 2
α, β, f, u = dados_de_entrada(exemplo)

# Define o número de elementos em que o domínio [0,1] é dividido
ne = 8

# Define a matriz LG
LG = monta_LG(ne)

# Define a dimensão do espaço aproximado Vₘ e o vetor EQ
m, EQ = monta_EQ(ne)

# Define a matriz EQoLG
EQoLG = EQ[LG]

# Exibe os valores ne, LG, m, EQ e EQoLG
display("ne = "); display(ne)
display("LG = "); display(LG)
display("m = "); display(m)
display("EQ = "); display(EQ)
display("EQoLG = "); display(EQoLG)

# Constrói a matriz K e o vetor F
K = monta_K(α, β, ne, m, EQoLG)
F = monta_F(f, ne, m, EQoLG)

# Exibe a matriz K (sistema linear) e o vetor F (vetor força)
display("Matriz K:")
display(K)
display("Vetor F:")
display(F)

# Resolve o sistema linear K * c = F para obter os coeficientes da solução aproximada u_h(x)
c = K \ F
display("Solução aproximada U:")
display(c)

# # Exibe a solução exata nos nos internos da malha
# display("Solução exata nos nós internos:")
# h = 1.0 / ne  # Comprimento de cada elemento finito
# display(u.(h:h:1-h))

# # # Calcula o erro entre a solução exata e a solução aproximada
# erro = erro_norma_L2(u, [c;0], ne, EQoLG)
# display("Erro L²:")
# display(erro)

# # Define a discretização da malha com N nós internos, incluindo os pontos de fronteira
# malha = 0:h:1

# # Gera o gráfico comparando a solução aproximada com a solução exata
# plt = plot(0:0.01:1, u.(0:0.01:1), label="Exata", lw=3, title="Comparação: Solução Exata vs. Aproximada")
# plot!(plt, malha, [0; c; 0], label="Aproximada", lw=3, linestyle=:solid, markershape=:circle) # O "!" adiciona ao gráfico existente
# xlabel!("x")  # Adiciona o rótulo ao eixo x

# # Exibe o gráfico final
# display(plt)