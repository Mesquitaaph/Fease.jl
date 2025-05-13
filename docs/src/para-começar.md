# Para começar

## Com projeto clonado (TESTAR EM OUTRO PC)
Verifique se o ambiente ativo é o `MyProject`. Para isso, no REPL da Julia, onde está escrito `julia>`, digite `]`. Se tudo der certo, aparecerá `(MyProject) pkg>` no lugar. Caso não apareça, execute o comando `activate MyProject` e deverá corrigir. Agora apague com o _backspace_.

```julia-repl
julia>

(xxxx) pkg> activate MyProject
  Activating new project at `~\MyProject`
```

Para verificar se o pacote está funcionando corretamente execute o comando `test` como apresentado abaixo.

```julia-repl
(MyProject) pkg> test MyProject
  Testing MyProject
  Status ...
```
Ao final deverá aparecer algo como
```julia-repl
    Testing Running tests...
Test Summary: | Pass  Total  Time
caso1D.jl     |    2      2  0.6s
Test Summary: | Pass  Total  Time
caso2D.jl     |    2      2  0.1s
    Testing MyProject tests passed 
```

A partir disso, é possível utilizar o MyProject incluindo o trecho abaixo no topo de seu código
```julia
using MyProject
```

```julia
example = 1
alpha, beta, gamma, a, b, u, u_x, f = examples(example)
run_values = RunValues(alpha, beta, gamma, f, u)
```

```julia
Nx1 = ne = 2^3
```

```julia
baseType = BaseTypes.linearLagrange
base = monta_base(baseType, ne)
```

```julia
malha = monta_malha_1D_uniforme(Nx1, base, a, b)
```

```julia
K = montaK_geral(run_values, malha)

F = montaF_geral(run_values, malha)
```

```julia
C = K\F
```

```julia
display(C)
```