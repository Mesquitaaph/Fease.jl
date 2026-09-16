function format_num(n)
  units = ["ns", "\$\\mu\$s", "ms", "s"]
  unit = ceil(Int, trunc(Int, log10(n) + 1) / 3)

  exp_div3 = trunc(Int, trunc(Int, log10(n)) / 3)

  num_sized = n / exp10(exp_div3 * 3)
  size_limit = 7
  size = length(string(num_sized)) > size_limit ? size_limit : length(string(num_sized))
  return string(num_sized)[1:size] * units[unit]
end

function measure_func(func, args)
  bench = @benchmark $func(($args)...)
  return (
    worst = maximum(bench.times),
    best = minimum(bench.times),
    mean = mean(bench.times),
    allocs = bench.allocs,
    memory = bench.memory
  )
end

function display_fieldnames(obj)
  return fieldnames(typeof(obj))
end

function test_revise() # Nos testes, verificar se essa saida é true
  return true
end

function showEQ(Nx1::Int64, Nx2::Int64, n_dir::Int64, EQ)
  EQs = []  
  for dir in 1:n_dir
    EQ_matrix = transpose(reshape(EQ[:, dir], (Nx1+1, Nx2+1)))
    append!(EQs, [EQ_matrix])
  end

  final_EQ = fill([], (Nx1+1, Nx2+1))
  for i in 1:(Nx1+1)*(Nx2+1)
    item = []
    for dir in 1:n_dir
      append!(item, EQs[dir][i])
    end

    final_EQ[i] = item
  end
  
  display(reverse(reverse(final_EQ), dims=2))
end

function getQuad(x, y, LG, elem)
  idx00 = LG[elem, 1]
  idx01 = LG[elem, 2]
  idx11 = LG[elem, 3]
  idx10 = LG[elem, 4]

  X = [x[idx00], x[idx10], x[idx11], x[idx01]]
  Y = [y[idx00], y[idx10], y[idx11], y[idx01]]

  return X, Y
end

function drawQuadInFig(fig, xVertices, yVertices)
  for i in 1:3
    xPair = [xVertices[i], xVertices[i+1]]
    yPair = [yVertices[i], yVertices[i+1]]
    linesegments!(fig[1, 1], xPair, yPair, color = :navy)
  end

  xPair = [xVertices[4], xVertices[1]]
  yPair = [yVertices[4], yVertices[1]]
  linesegments!(fig[1, 1], xPair, yPair, color = :navy)

  return fig
end

function drawGrid(xValues, yValues, numX, numY, LG, windowX, windowY)
  LGᵗ = transpose(LG)

  f = Figure(size = (500, 500))
  Axis(f[1, 1], limits = (0, windowX, 0, windowY))

  for i in 1:numY
    for j in 1:numX
      elem = (i - 1) * (numX) + j

      xVertices, yVertices = getQuad(xValues, yValues, LGᵗ, elem)
      drawQuadInFig(f, xVertices, yVertices)
    end
  end

  return f
end
export drawGrid

function drawGridFromFile(elemX, elemY, coordsFileName, LGFileName, xWindowLim, yWindowLim)
  malhaFile = readdlm(coordsFileName)

  X, Y = malhaFile[:, 2], malhaFile[:, 3]

  LG = readdlm(LGFileName)[:, 1:4]
  LG = convert(Matrix{Int}, LG)
  LG = LG .+ 1 # Indices 1-based (julia)

  return drawGrid(X, Y, elemX, elemY, transpose(LG), xWindowLim, yWindowLim)
end
export drawGridFromFile

function calcula_erro(passo, F_def, s1, s2, τ, pts_fronteira, ne)
  max_F11 = 0.0
  min_F11 = 2.0
  sum_F11 = 0.0
  sum_det = 0.0

  for F in F_def
    sum_det += det(F)
    sum_F11 += F[1, 1]
    max_F11 = (F[1, 1] > max_F11) ? F[1, 1] : max_F11
    min_F11 = (F[1, 1] < min_F11) ? F[1, 1] : min_F11
  end

  λ₁_analitico = (1.0 - τ^2 * (1.0 / (s1 - s2)^2))^(-1.0 / 4)
  k = sqrt(λ₁_analitico^4 - 1.0)

  p1 = pts_fronteira[1]
  p2 = pts_fronteira[2]
  p4 = pts_fronteira[4]

  OA_num = p2[1] - p1[1]
  OB_num = sqrt((p4[1]^2 - p1[1]^2) + (p4[2]^2 - p1[2]^2))

  OA_analitico = λ₁_analitico^2
  #F_11 médio
  λ₁_numerico = sum_F11 / ne
  det_medio = sum_det / ne

  erro_relativo_percentual = 100 * abs((λ₁_numerico - λ₁_analitico) / λ₁_analitico)
  razao_segmentos = 100 * abs((OA_num - OB_num) / OA_num)
  infos_principais = [passo,
    λ₁_analitico,
    λ₁_numerico,
    min_F11,
    max_F11,
    det_medio,
    erro_relativo_percentual,
    razao_segmentos]::Vector{Float64}

  infos_secundarias = [k, OA_analitico]

  return infos_principais, infos_secundarias
end
export calcula_erro