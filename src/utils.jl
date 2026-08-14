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