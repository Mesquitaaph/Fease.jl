function plot_solucao_aproximada(c::Array, malha::Malha, gif::Bool = false)
  uh = monta_u_aproximada(c, malha)

  return plot_solucao_aproximada(uh, malha, gif)
end

function plot_solucao_aproximada(uh::Function, malha::Malha, gif::Bool = false)
  (; a, b, n_dim) = malha

  X = ()
  for dim in 1:n_dim
    X = (X..., range(a[dim], b[dim], length = 40))
  end

  ST = (:path, :surface)

  # create a plot with 3 subplots and a custom layout
  p = plot(X..., uh, st = ST[n_dim])
  if !gif
    return display(p)
  end

  n = 100
  @gif for i in range(0, stop = 360, length = n)

    # induce a slight oscillating camera angle sweep, in degrees (azimuth, altitude)
    plot!(p[1], camera = (i, 40))
    display(p)
  end
end

function plot_malha_2D(
    malha::Malha
)
  # Cria uma nova figura
  fig = Plots.plot(legend = false, aspect_ratio = :equal,
    xticks = 0:0.25:1,  # Defina os ticks do eixo x
    yticks = 0:0.25:1)  # Defina os ticks do eixo y

  (; coords, LG, ne) = malha
  (X₁, X₂) = coords
  # Adiciona os nós da malha
  Plots.scatter!(X₁, X₂, markersize = 4, color = :blue)

  # Adiciona as linhas de contorno dos elementos
  for e in 1:ne
    # Obtém os índices dos nós do elemento `e`
    i₁, i₂, i₃, i₄ = LG[:, e]

    # Adiciona as linhas de contorno do elemento `e`
    Plots.plot!([X₁[i₁], X₁[i₂], X₁[i₄], X₁[i₃], X₁[i₁]],
      [X₂[i₁], X₂[i₂], X₂[i₄], X₂[i₃], X₂[i₁]], color = :black)
  end

  return display(fig)
end
