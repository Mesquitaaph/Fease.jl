# Para começar a desenvolver
Neste capítulo esclareceremos os passos necessários para começar a desenvolver e contribuir com o pacote.

## Clonando o projeto
Em um terminal com o diretório de trabalho em uma pasta de sua escolha, execute os comandos abaixo (Windows ou Linux)
```shell
git clone https://github.com/Mesquitaaph/MyProject.jl.git
cd MyProject
```
Em seguida, você pode abrir o projeto com o VSCode com mais um comando
```shell
code .
```

## Durante o desenvolvimento... (Essa parte talvez entre nas diretrizes de dev colab)
Acredito que a melhor forma de testar uma nova implementação seria incluindo um novo arquivo com exemplo de como utilizá-la, na pasta `exemplos`, na raíz do projeto. Inicialmente existem dois arquivos: `caso1D.jl` e `caso2D.jl`, cada um com um exemplo distinto.

## Adicionando o pacote
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

## Instalando e configurando o Revise.jl