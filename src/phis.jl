function ϕ_geral(P...)
  # Acho que vou alterar essa função pra funcionar como a ∇ϕ_geral.
  # Ela tem um funcionamento complicado, apesar de parecido com o cálculo das combinações dos pontos de Gauss
  # Também não sei se funciona para outras bases diferentes da lagrange linear
  # Pretendo trocar ela para funcionar como a ϕ_1D abaixo
  n_dim = length(P)
  phis = 2.0^(-n_dim) * ones(Float64, 2^n_dim)
  for d in 1:n_dim
    for n_phi in 1:(2^n_dim)
      sinal = floor((2.0^(d - 1) + n_phi - 1) * 2.0^-(d - 1))
      phis[n_phi] *= 1 + ((-1.0)^sinal) * P[d]
    end
  end

  return (; phis)
end

function ϕ_1D(P...)
  # Ela é preparada, não intencionalmente, para receber números de dimensão igual a 1 ou a 2
  n_dim = length(P)
  if n_dim == 1
    phis = 2.0^-n_dim * [
      1 - P[1],
      1 + P[1]
    ]
    return (; phis)
  elseif n_dim == 2
    return 2.0^-n_dim * [
      (1 - P[1]) * (1 - P[2]),
      (1 + P[1]) * (1 - P[2]),
      (1 - P[1]) * (1 + P[2]),
      (1 + P[1]) * (1 + P[2])
    ]
  else
    return Nothing
  end
end

function ϕ_2D(ξ₁::Float64, ξ₂::Float64)::Vector{Float64}
  # É apenas a definição das ϕ lagrange linear para o caso 2D
  return [
    (1 - ξ₁) * (1 - ξ₂) / 4,
    (1 + ξ₁) * (1 - ξ₂) / 4,
    (1 - ξ₁) * (1 + ξ₂) / 4,
    (1 + ξ₁) * (1 + ξ₂) / 4
  ]
end

function ∇ϕ_1D(P...)
  # É apenas a definição das dϕ lagrange linear para o caso 1D
  # Teste de número de dimensão (n_dim) é resquício e deverá ser apagado
  n_dim = length(P)
  if n_dim == 1
    dphis = 2.0^-n_dim * [
      -1,
      1
    ]
    return (; dphis)
  else
    return error("Dimensão não implementada")
  end
end

function ∂ϕ_∂ξ₁(ξ₂::Float64)::Vector{Float64}
  # É apenas a definição das ∂ϕ_∂ξ₁ lagrange linear para o caso 2D
  return [-(1 - ξ₂) / 4, (1 - ξ₂) / 4, -(1 + ξ₂) / 4, (1 + ξ₂) / 4]
end

function ∂ϕ_∂ξ₂(ξ₁::Float64)::Vector{Float64}
  # É apenas a definição das ∂ϕ_∂ξ₂ lagrange linear para o caso 2D
  return [-(1 - ξ₁) / 4, -(1 + ξ₁) / 4, (1 - ξ₁) / 4, (1 + ξ₁) / 4]
end

function ∇ϕ_2D(ξ₁::Float64, ξ₂::Float64)
  # Retorna as duas como um gradiente, de fato
  return (∂ϕ_∂ξ₁(ξ₂), ∂ϕ_∂ξ₂(ξ₁))
end

function ∇ϕ_geral(P...)
  # Apenas retorna ∇ϕ respectivo ao número da dimensão
  n_dim = length(P)

  if n_dim == 1
    return ∇ϕ_1D(P...)
  elseif n_dim == 2
    return ∇ϕ_2D(P...)
  else
    return error("Dimensão não implementada")
  end
end
