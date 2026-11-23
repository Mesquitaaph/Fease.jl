"""
    ϕ_geral(P...)

Calcula os valores das funções interpoladoras `ϕ` no ponto `P`.

# Parâmetros
- `P...`: Ponto onde as funções interpoladoras serão calculadas.

# Retorno
- `ϕ::Vector{Float64}`: Vetor com os valores das funções interpoladoras `ϕ` no ponto `P`.

# Exemplo
```@example

```
"""
function ϕ_geral(P...)
    n_dim = length(P)
    if n_dim == 1
        return ϕ_1D(P...)
    elseif n_dim == 2
        return ϕ_2D(P...)
    else
        return error("Dimensão não implementada")
    end
end

"""
    ϕ_1D(ξ)
    
Calcula os valores das funções interpoladoras `ϕ` no ponto `ξ`, unidimensional.

# Parâmetros
- `ξ`: Ponto onde as funções interpoladoras serão calculadas.

# Retorno
- `ϕ::Vector{Float64}`: Vetor com os valores das funções interpoladoras `ϕ` no ponto `ξ`.

# Exemplo
```@example

```
"""
function ϕ_1D(ξ)::Vector{Float64}
    # Retorna diretamente o Vector para o caso 1D
    return [
        (1 - ξ) / 2.0,
        (1 + ξ) / 2.0
    ]
end

"""
    ϕ_2D(ξ₁::Float64, ξ₂::Float64)::Vector{Float64}
    
Calcula os valores das funções interpoladoras `ϕ` no ponto `(ξ₁, ξ₂)`, bidimensional.

# Parâmetros
- `ξ₁::Float64`: Valor no primeiro eixo do ponto onde as funções interpoladoras serão calculadas.
- `ξ₂::Float64`: Valor no segundo eixo do ponto onde as funções interpoladoras serão calculadas.

# Retorno
- `ϕ::Vector{Float64}`: Vetor com os valores das funções interpoladoras `ϕ` no ponto `(ξ₁, ξ₂)`.

# Exemplo
```@example

```
"""
function ϕ_2D(ξ₁::Float64, ξ₂::Float64)::Vector{Float64}
    # É apenas a definição das ϕ lagrange linear para o caso 2D
    # Ordem anti-horária
    return [
        (1 - ξ₁) * (1 - ξ₂) / 4,
        (1 + ξ₁) * (1 - ξ₂) / 4,
        (1 + ξ₁) * (1 + ξ₂) / 4,
        (1 - ξ₁) * (1 + ξ₂) / 4]
end

export ϕ_2D

"""
    ∇ϕ_1D(ξ)

Calcula os valores das derivadas das funções interpoladoras `ϕ` no ponto `ξ`.

# Parâmetros
- `ξ...`: Ponto onde as funções interpoladoras serão calculadas.

# Retorno
- `∇ϕ::Vector{Float64}`: Vetor com os valores das derivadas das funções interpoladoras `ϕ` no ponto `ξ`.

# Exemplo
```@example

```
"""
function ∇ϕ_1D(ξ)
    # É apenas a definição das dϕ lagrange linear para o caso 1D
    dphis = [
        -1/2.0,
         1/2.0
    ]
    return (; dphis)
end
export ∇ϕ_1D

"""
    ∂ϕ_∂ξ₁(ξ₂::Float64)::Vector{Float64}

Definição das ∂ϕ_∂ξ₁ lagrange linear para o caso 2D

# Parâmetros
- `ξ₂::Float64`: Valor no segundo eixo do ponto onde as funções interpoladoras serão calculadas.

# Retorno
- `∂ϕ_∂ξ₁::Vector{Float64}`: Vetor com os valores das derivadas das funções interpoladoras `ϕ` no ponto `P`, no primeiro eixo.

# Exemplo
```@example

```
"""
function ∂ϕ_∂ξ₁(ξ₂::Float64)::Vector{Float64}
    # É apenas a definição das ∂ϕ_∂ξ₁ lagrange linear para o caso 2D
    # Ordem anti-horária
    return [-(1 - ξ₂) / 4, (1 - ξ₂) / 4, (1 + ξ₂) / 4, -(1 + ξ₂) / 4]
end

"""
    ∂ϕ_∂ξ₂(ξ₁::Float64)::Vector{Float64}

Definição das ∂ϕ_∂ξ₂ lagrange linear para o caso 2D

# Parâmetros
- `ξ₁::Float64`: Valor no primeiro eixo do ponto onde as funções interpoladoras serão calculadas.

# Retorno
- `∂ϕ_∂ξ₂::Vector{Float64}`: Vetor com os valores das derivadas das funções interpoladoras `ϕ` no ponto `P`, no segundo eixo.

# Exemplo
```@example

```
"""
function ∂ϕ_∂ξ₂(ξ₁::Float64)::Vector{Float64}
    # É apenas a definição das ∂ϕ_∂ξ₂ lagrange linear para o caso 2D
    # Ordem anti-horária
    return [-(1 - ξ₁) / 4, -(1 + ξ₁) / 4, (1 + ξ₁) / 4, (1 - ξ₁) / 4]
end

"""
    ∇ϕ_2D(ξ₁::Float64, ξ₂::Float64)

Calcula os valores do gradiente das funções interpoladoras `ϕ` no ponto `(ξ₁, ξ₂)`.

# Parâmetros
- `ξ₁::Float64`: Valor no primeiro eixo do ponto onde as funções interpoladoras serão calculadas.
- `ξ₂::Float64`: Valor no segundo eixo do ponto onde as funções interpoladoras serão calculadas.

# Retorno
- `∇ϕ::Tuple`: Tupla com os valores do gradiente das funções interpoladoras `ϕ` no ponto `(ξ₁, ξ₂)`.

# Exemplo
```@example

```
"""
function ∇ϕ_2D(ξ₁::Float64, ξ₂::Float64)
    # Retorna as duas como um gradiente, de fato
    return (∂ϕ_∂ξ₁(ξ₂), ∂ϕ_∂ξ₂(ξ₁))
end
export ∇ϕ_2D

"""
    ∇ϕ_geral(P...)

Calcula os valores do gradiente das funções interpoladoras `ϕ` no ponto `P`.

# Parâmetros
- `P...`: Ponto onde as funções interpoladoras serão calculadas.

# Retorno
- `∇ϕ::Tuple`: Tupla com os valores do gradiente das funções interpoladoras `ϕ` no ponto `P`.

# Exemplo
```@example

```
"""
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
export ∇ϕ_geral
