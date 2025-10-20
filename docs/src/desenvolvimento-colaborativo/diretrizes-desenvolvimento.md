# Diretrizes para o Desenvolvimento Colaborativo

Aqui apresentaremos o que deve ser seguido durante a implementação de novas funcionalidades, melhorias ou correções com o objetivo de manter o projeto sustentável e escalável, englobando desenvolvedores com experiências de programação em diferentes níveis.

# Durante o desenvolvimento...

A melhor forma de testar uma nova implementação é incluindo um novo arquivo na pasta src/testes, como apresentado em [Para Começar a Desenvolver](../para-comecar/para-começar-a-desenvolver.md).

Uma vez finalizada, se sua implementação incluir algo novo, que seja necessário um tutorial, como por exemplo, um problema que não existe, deve-se incluir um exemplo base na pasta example na raíz do projeto. Inicialmente existem dois arquivos: caso1D.jl e caso2D.jl, cada um com um exemplo distinto.

# Legibilidade

A legibilidade do código é definida por alguns pontos:

## Idioma do código

Apesar da documentação ser escrita em português, nosso código deve ser escrito em **inglês**.

## Práticas de Escrita - Convenção de Nomenclatura

Deve-se utilizar `snake_case` para nome de funções e variáveis, por exemplo, `test_revise` e `CamelCase` para tipos implementados e arquivos, por exemplo, `EquationTerms`.

## Nomes de funções, variáveis, tipos e nomes de arquivo

Funções e variáveis devem ser nomeadas de forma clara e descritiva, facilitando a leitura e compreensão de quem vai trabalhar com sua implementação. Sempre que possível use palavras inteiras. Por exemplo: escrevi uma função que monta e retorna uma malha 2D. Nomearei-a como `assemble_mesh_2D`.

# Identação, espaçamento e outras definições de escrita

Essas outras escolhas de estilo de escrita são aplicadas automaticamente pelo _formatter_ configurado no projeto. No VSCode, todos os arquivos são auto-formatados ao salvar, não sendo necessário se preocupar com isso.

Caso não esteja utilizando o VSCode, deve-se executar

```julia-repl
julia> using Fease
julia> format(".")
```

para formatar todos os arquivos modificados.
