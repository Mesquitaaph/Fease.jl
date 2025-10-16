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

function convergence_test!(E, NE, run_values, baseType, a, b, n_dims)
  fill!(E, 0.0)
  dE = similar(E)
  (; u) = run_values

  for i in 1:lastindex(NE)
    NX = []
    ne = 1
    for dim in 1:n_dims
      e = NE[i]
      push!(NX, e)
      ne *= e
    end

    base = monta_base(baseType, ne)

    malha = Nothing

    if n_dims == 1
      malha = monta_malha_1D_uniforme(base, NX[1], a, b)
    elseif n_dims == 2
      malha = monta_malha_2D_uniforme(base, NX..., a, b)
    end

    C = solve_sys_poisson(run_values, malha)

    E[i], dE[i] = erro_L2(malha, u, C)
  end
end
