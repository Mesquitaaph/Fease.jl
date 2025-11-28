# Para começar a utilizar

Neste capítulo esclareceremos os passos necessários para começar a utilizar o pacote.

## Instalando o pacote

O pacote pode ser instalado em seu ambiente de duas formas distintas: através dos gerenciador de pacotes de Julia ou através do projeto clonado. As apresentaremos a seguir.

#### Através do gerenciador de pacotes de Julia

Ainda não subimos no gerenciador de pacotes.

#### Através do projeto clonado

Ao iniciar o REPL da Julia, você verá escrito `julia>`. Para utilizar o Fease.jl é preciso adicioná-lo ao seu ambiente. Para isso, digite `]` na linha de comando para entrar no modo _package_. Deverá aparecer algo como `(nome_do_ambiente) pkg>` no lugar. Em seguida, adicione-o como um pacote de desenvolvimento digitando

```julia-repl
(nome_do_ambiente) pkg> dev C:\caminho\para\o\Fease
  Resolving package versions...
```

Para verificar se o pacote está funcionando corretamente execute o comando `test` como apresentado abaixo.

```julia-repl
(nome_do_ambiente) pkg> test Fease
  Testing Fease
  Status ...
```

Ao final deverá aparecer algo como

```julia-repl
    Testing Running tests...
Test Summary: | Pass  Total  Time
caso1D.jl     |    2      2  0.6s
Test Summary: | Pass  Total  Time
caso2D.jl     |    2      2  0.1s
    Testing Fease tests passed
```

## Utilizando o pacote

A partir disso, é possível utilizar o Fease incluindo o trecho abaixo no topo de seu código

```julia
using Fease
```

O fluxo de implementação inicia definindo o conjunto de funções base. Primeiro declaramos a quantidade de elementos finitos. Em seguida, o tipo das funções base.

```julia
Nx1 = 2^3

baseType = BaseTypes.linearLagrange
```

Com a base definida, construímos a malha utilizando algumas funções já implementadas, como `monta_malha_1D_uniforme`. Para esta, precisamos declarar os pontos inicial e final do intervalo $[a, b]$.

```julia
# Define extremos do intervalo da malha.
a, b = 0, 1

# Define a malha com os valores atribuídos acima.
malha = monta_malha_1D_uniforme(baseType, Nx1, a, b)
```

O próximo passo é definir os valores do problema. Para este tutorial, podemos utilizar valores de problemas exemplo do Caso 1D

```julia
# Define os parâmetros da equação a ser resolvida.
example = 1
run_values = examples_1D(example)
(; α, β, f) = run_values
```

Em seguida, montamos uma função que referencia o operador bilinear $a(u,v)$ com parâmetros obtidos de `run_values`.

```julia
# Define o pseudo operador linear a(u,v).
function pseudo_a(termos_equacao)
  (; ∇u, ∇v, u, v) = termos_equacao

  return β * dot(u, v) + α * dot(∇u, ∇v)
end
```

Enfim, com a `f` e a malha definidas e o operador $a(u, v)$ referenciado, basta resolver o sistema, digitando

```julia
# Monta e resolve o sistema linear relacionado a esta equação.
c = solve_sys(f, malha, pseudo_a)
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

Para mais informações de uso, vá para [Tutoriais](../tutoriais/index.md)
