# Implementação Base

Aqui esclareceremos como foi pensada a implementação base deste projeto.

## Generalização

## Estruturas

No código existem estruturas e tipos que são essenciais para seu funcionamento e entendimento.

### Mesh

Aqui se encontram valores que definem uma malha. Nela se encontram, por exemplo, o número de subdivisões para cada eixo, número de elementos totais, as coordenadas de cada um dos nós, sua dimensão, a matriz LG, o vetor EQ, entre outros.

Existem funções para montagens específicas de malha. Um exemplo é a `assemble_uniform_mesh_1D` que retorna uma variável do tipo `Mesh` para um intervalo uniforme discreto dentro de $[a, b],\ a<b$.

Manter todos esses valores agrupados, torna-os facilmente acessíveis.

### RunValues

Essa estrutura é utilizada para agrupar valores normalmente utilizados nos problemas, ou seja, parâmetros, coeficientes, a $f$ e a solução $u$ esperada. É o que retorna quando se chama funções "exemplos", localizadas em `examples.jl`. Seu uso é opcional, recomendado para montar esses exemplos prontos.

### EquationTerms

Para determinar os termos a serem utilizados no operador bilinear $a(u,v)$, é necessário o `EquationTerms`. É uma forma de indicar como $u$, $v$, $\nabla u$, $\nabla v$... são aplicados em $a(u, v)$.

## Montando operador $a(u, v)$ e aplicando-o na montagem do sistema linear

O operador pode se encontrar, por exemplo, na seguinte forma:

$$a(u, v) = \beta \int{ u(x)\cdot v(x)d\Omega}  + \alpha \int{ \nabla{u}(x)\cdot \nabla{v}(x) d\Omega}$$

No entanto, o que utilizamos é uma forma sem a integração, ou seja, um pseudo operador, por exemplo, `pseudo_a(u,v) = β * dot(u, v) + α * dot(∇u, ∇v)`.
Nesse caso, α e β são constantes escolhidas arbitrariamente, mas `u`, `v`, `∇u` e `∇v` são referenciadas diretamente de `EquationTerms`.

Então, a função que utilizamos para definir esse pseudo operador, na verdade, é:

```julia
function pseudo_a(equation_terms::EquationTerms)
  (; ∇u, ∇v, u, v) = equation_terms

  return β * dot(u, v) + α * dot(∇u, ∇v)
end
```

O argumento `equation_terms` contém todos os termos possíveis de serem utilizados e pegamos apenas os que queremos em `(; ∇u, ∇v, u, v) = equation_terms`. Essa nomenclatura `(; var1, var2, var3,...)` no lado esquerdo da atribuição permite extrair os campos de `EquationTerms` através de seus nomes, sem precisar de uma ordem específica. Se quisermos um outro, por exemplo, `x`, podemos digitar `(; ∇u, ∇v, u, x, v) = equation_terms`.

Já o retorno da função é o que será aplicado na montagem das `K locais` do problema aproximado. Temos que uma entrada da matriz $K^e$ local é

$$K^{e}_{a,b} = \int \alpha \nabla{\varphi^{e}_{b}(x(\xi))} \cdot \nabla{\varphi^{e}_{a}(x(\xi))} + \beta \varphi^{e}_{b}(x(\xi)) \cdot \varphi^{e}_{a}(x(\xi))\ |J(\xi)|d\xi_{1}d\xi_{2}$$

Essa integração é feita com a Quadratura Gaussiana, assim, para cada ponto de gauss $i$, temos a expressão

$$w_i(\alpha \nabla{\varphi^{e}_{b}(x(\xi))} \cdot \nabla{\varphi^{e}_{a}(x(\xi))} + \beta \varphi^{e}_{b}(x(\xi)) \cdot \varphi^{e}_{a}(x(\xi)))|J(\xi)|$$

Portanto, se substituirmos $\nabla{\varphi^{e}_{b}(x(\xi))}$ por `∇u`, $\nabla{\varphi^{e}_{a}(x(\xi))}$ por `∇v`, $\varphi^{e}_{b}(x(\xi))$ por `u` e $\varphi^{e}_{b}(x(\xi))$ por `v`, temos `wᵢ*(α*∇u⋅∇v + β*u⋅v)|J(ξ)|`. Ou seja, a expressão entre parênteses é examente o que a função `pseudo_a` retorna.

Isso na montagem da `K local` é feito com o trecho de código a seguir

```julia
  equation_terms = EquationTerms(
    ϕᵉ_b,
    ϕᵉ_a,
    ∇ϕᵉ_b,
    ∇ϕᵉ_a
  )
  sum = pseudo_a(equation_terms)
  "sum = β * dot(ϕᵉ_b, ϕᵉ_a) + α * dot(∇ϕᵉ_b, ∇ϕᵉ_a)"

  Kᵉ[a, b] += wᵢ * sum * detJ
```

Assim, é possível definir uma referência ao operador $a(u, v)$, com os termos apresentados, sem precisar mexer no código da montagem da `K local`.

### Adicionando um novo termo

Adicionar um novo termo vem da ideia de tornar flexível a implementação de novos problemas. Necessariamente tem que ser calculado e obtido dentro da montagem da `K local`, assim como `∇ϕᵉ_b`, por exemplo. Um exemplo pronto de como incluir um termo novo se encontra no [tutorial](../tutoriais/novo-termo.md)
