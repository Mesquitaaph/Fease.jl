function erro_L2(malha::Malha, u, C)
  npg = 5

  (; ne, neq, n_dim, base) = malha

  # Avalia as ϕ e ∇ϕ nos pontos de Gauss
  ϕξ, P, W = quadratura_ϕ(base, npg, n_dim)
  ∇ϕξ, P, W = quadratura_∇ϕ(base, npg, n_dim)

  d = [C..., 0]

  EL2 = 0.0
  for e in 1:ne
    eqs_idx, Xᵉ = elem_coords(malha::Malha, e::Int)

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

      approx = 0
      for a in 1:(2^n_dim)
        i = eqs_idx[a]
        approx += d[i] * ϕᵉ[a]
      end

      @inbounds EL2 += WW * (u(x...) - approx)^2 * detJ
    end
  end

  EL2 = sqrt(EL2)

  EH01 = 0.0
  # EH01 = sqrt(h/2 * sum(W' * ((u_x.(xPTne) - (2/h .* dphiP * cEQoLG)).^2)))

  return EL2, EH01
end

"""
    convergence_test!(E, NE, n_dims, monta_malha, pseudo_a, f, u)

Realiza o estudo da convergência do erro para a solução aproximada da equação com solução exata `u`.

# Parâmetros
- `E::Vetor{Float64}`: O vetor onde serão inseridos os erros caculados para cada um das quantidades de elementos finitos.
- `NE::Vetor{Float64}`: O vetor que contém as cada uma das quantidades de elementos finitos.
- `n_dim::Int`: Número de dimensões do espaço.
- `monta_malha::Function`: Uma função que recebe os números de sub-intervalos da malha e a retorna montada.
- `pseudo_a`: A referência ao operador bilinear a(u,v).
- `f::Function`: Função que representa o lado direito da equação do problema.
- `u::Function`: A solução analítica da equação.

# Retorno
Altera o vetor `E`.

# Exemplo
```@example
  # Define o tipo da base de funções interpoladoras do subespaço aproximado Vₘ.
  baseType = BaseTypes.linearLagrange

  # Define extremos do intervalo da malha.
  a = (0.0, 0.0)
  b = (1.0, 1.0)

  # Define os parâmetros da equação a ser resolvida.
  α = 1.0
  β = 1.0
  f = (x₁, x₂) -> (2 * α * π^2 + β) * sin(π * x₁) * sin(π * x₂)
  u = (x₁, x₂) -> sin(π * x₁) * sin(π * x₂)

  # Define o pseudo operador linear a(u,v).
  function pseudo_a(termos_equacao)
    (; ∇u, ∇v, u, v) = termos_equacao

    return β * dot(u, v) + α * dot(∇u, ∇v)
  end

  # Define a função que constrói uma malha que depende só do número de sub-intervalos nos dois eixos.
  function monta_malha(NX)
    return monta_malha_2D_uniforme(baseType, NX..., a, b)
  end

  # Estudo será realizado com número de elementos em potências de 2.
  errsize = 7
  NE = 2 .^ [2:1:errsize;]
  H = 1 ./ NE
  E = zeros(length(NE))

  # Calcula o erro para cada uma das quantidades de elementos finitos.
  convergence_test!(E, NE, 2, monta_malha, pseudo_a, f, u)
```
"""
function convergence_test!(E, NE, n_dim, monta_malha, pseudo_a, f, u)
  fill!(E, 0.0)
  dE = similar(E)

  for i in 1:lastindex(NE)
    NX = []
    ne = 1
    for dim in 1:n_dim
      e = NE[i]
      push!(NX, e)
      ne *= e
    end

    malha = monta_malha(NX)

    C = solve_sys(f, malha, pseudo_a)

    E[i], dE[i] = erro_L2(malha, u, C)
  end
end
