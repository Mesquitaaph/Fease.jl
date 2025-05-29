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

## Desenvolvendo com VSCode
Recomendamos que seja instalada a extensão de Julia em sua instalação do VSCode. Com o auxílio deste 
[passo a passo](https://code.visualstudio.com/docs/getstarted/extensions)
, pesquise por `julialang.language-julia` e clique em "instalar".

Agora você deve iniciar o REPL de Julia. Com um arquivo `.jl` do projeto, pode fazer isso apertando as teclas `Ctrl+Shift+P` e digitar `Start REPL`.

Ao iniciar o REPL da Julia, você verá escrito `julia>` no terminal integrado do VSCode. Para testar o MyProject.jl digite `]` na linha de comando para entrar no modo _package_. Deverá aparecer `(MyProject) pkg>` no lugar. Para verificar se o pacote está funcionando corretamente execute o comando `test` como apresentado abaixo.

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

## Instalando e configurando o Revise.jl
O Revise.jl é um pacote a ser instalado globamente que auxilia no desenvolvimento de pacotes Julia. Ele atualiza o pacote local automaticamente ao salvar o arquivo.
Isso reduz a quantidadade de vezes que precisamos reiniciar o ambiente/REPL.

Ainda com o terminal integrado do VSCode aberto, entre no modo _package_ e digite `activate`:
```julia-repl
(MyProject) pkg> activate
```
Deverá mudar de `(MyProject) pkg>` para `(@v1.11) pkg>` (ou para a versão de Julia que estiver instalada em sua máquina). Com isso, prossiga:
```julia-repl
(@v1.11) pkg> add Revise
```
Adicionando-o com sucesso, o Revise.jl estará instalado globalmente.

Para configurá-lo, abra um outro terminal do seu sistema operacional e execute:
- Terminal do Linux
```shell
mkdir -p ~/.julia/config/ && echo "using Revise" >> ~/.julia/config/startup.jl
```
- Cmd do Windows
```shell
mkdir %userprofile%\.julia\config && echo using Revise >> %userprofile%\.julia\config\startup.jl
```

Agora reabra o VSCode e inicie o REPL. Execute:
```julia-repl
julia> using MyProject
```

## Testando mudanças no pacote
Na raíz do projeto, dentro da pasta `src/testes`, crie um arquivo Julia, por exemplo `teste.jl` e inclua-o no arquivo `include_testes.jl` digitando, ao final do arquivo: 

```julia
include("teste.jl")
```

De volta ao arquivo criado, digite:
```julia
test_revise()
```

Para executar essa função clique no botão de _play_ no canto superior direito da janela do VSCode para executar o arquivo inteiro ou clicando na linha com o trecho do código e apertando no teclado `Ctrl+Enter`. O resultado será apresentado no terminal, rodando a função com a primeira opção, ou na linha, com a segunda opção. Deverá aparecer `true`.

Agora segure a tecla `Ctrl` e clique na função `test_revise`. Este atalho abrirá o arquivo utils.jl, onde encontra sua definição, na linha certa. Troque seu retorno de `true` para `false`:
```julia
function test_revise()
  return false
end
```

Volte ao arquivo `teste.jl` e execute novamente a função. Se desta vez aparecer `false`, o `Revise` está funcionando e, assim, seu ambiente estará configurado e pronto para você desenvolver!