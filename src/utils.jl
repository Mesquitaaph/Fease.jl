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
