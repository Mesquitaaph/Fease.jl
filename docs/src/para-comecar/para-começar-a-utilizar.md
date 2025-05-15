# Para começar a utilizar
Neste capítulo esclareceremos os passos necessários para começar a utilizar o pacote.

## Instalando o pacote
O pacote pode ser instalado em seu ambiente de duas formas distintas: através dos gerenciador de pacotes de Julia ou através do projeto clonado. As apresentaremos a seguir.

#### Através do gerenciador de pacotes de Julia
Ainda não subimos no gerenciador de pacotes.

#### Através do projeto clonado
Ao iniciar o REPL da Julia, você verá escrito `julia>`. Para utilizar o MyProject.jl é preciso adicioná-lo ao seu ambiente. Para isso, digite `]` na linha de comando para entrar no modo _package_. Deverá aparecer algo como `(nome_do_ambiente) pkg>` no lugar. Em seguida, adicione-o como um pacote de desenvolvimento digitando

```julia-repl
(nome_do_ambiente) pkg> dev C:\caminho\para\o\MyProject
  Resolving package versions...
```

Para verificar se o pacote está funcionando corretamente execute o comando `test` como apresentado abaixo.

```julia-repl
(nome_do_ambiente) pkg> test MyProject
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

## Utilizando o pacote
A partir disso, é possível utilizar o MyProject incluindo o trecho abaixo no topo de seu código
```julia
using MyProject
```

O fluxo de implementação inicia definindo o conjunto de funções base. Primeiro declaramos a quantidade de elementos finitos. Em seguida, o tipo das funções base. Enfim, chamamos o método `monta_base` para a obtermos nossa base.

```julia
ne = 2^3

baseType = BaseTypes.linearLagrange
base = monta_base(baseType, ne)
```

Com a base definida, construímos a malha utilizando algumas funções já implementadas, como `monta_malha_1D_uniforme`. Para esta, precisamos declarar os pontos inicial e final do intervalo $[a, b]$.
```julia
a = 0
b = 1
malha = monta_malha_1D_uniforme(ne, base, a, b)
```

O próximo passo é definir os valores do problema. Para este tutorial, podemos utilizar valores de problemas exemplo do Caso 1D
```julia
example = 1
(; alpha, beta, gamma, u, f) = examples(example)
run_values = RunValues(alpha, beta, gamma, f, u)
```

Em seguida, montamos a matriz $K$ e o vetor $F$ do sistema linear $Kc = F$. As duas formas de montagem são praticamente iguais, apenas mudando a função a ser chamada.
```julia
K = montaK_geral(run_values, malha)

F = montaF_geral(run_values, malha)
```

Enfim, com a $K$ e $F$ definidas, podemos utilizar o operador padrão de resolução de sistema linear `\` de Julia.
```julia
c = K\F
```

Obtendo o resultado
```julia
7-element Vector{Float64}:
 0.11108398583873315
 0.21339432867780742
 0.29760326069653276
 0.3532397670750393
 0.36802825874249573
 0.32711327201474155
 0.2121212287826553
```