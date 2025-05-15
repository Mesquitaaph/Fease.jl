using Documenter, MyProject

DocMeta.setdocmeta!(MyProject, :DocTestSetup, :(using MyProject); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers

function nice_name(file)
  file = replace(file, r"^[0-9]*-" => "")
  if haskey(page_rename, file)
    return page_rename[file]
  end
  return splitext(file)[1] |> x -> replace(x, "-" => " ") |> titlecase
end

Para_Comecar = "Para Começar" => [
  "A Utilizar" => "para-comecar/para-começar-a-utilizar.md",
  "A Desenvolver" => "para-comecar/para-começar-a-desenvolver.md",
]

Dev_Colab = "Desenvolvimento Colaborativo" => [
  "Diretrizes de Desenvolvimento" => "desenvolvimento-colaborativo/diretrizes-desenvolvimento.md",
  "Implementação Base" => "desenvolvimento-colaborativo/implementacao-base.md",
  "Ferramentas" => "desenvolvimento-colaborativo/ferramentas.md",
]

makedocs(;
  modules = [MyProject],
  doctest = true,
  linkcheck = true, # Rely on Lint.yml/lychee for the links
  authors = "Raphael Mesquita <raphaelfcm@ic.ufrj.br> and contributors",
  repo = "https://github.com/Mesquitaaph/MyProject.jl/blob/{commit}{path}#{line}",
  sitename = "MyProject.jl",
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical = "https://Mesquitaaph.github.io/MyProject.jl",
    # assets = ["assets/style.css"],
  ),
  pages = [
    "Início" => "index.md",
    Para_Comecar,
    "Método de Elementos Finitos" => "metodo-elementos-finitos.md",
    "Tutoriais" => "tutoriais.md",
    Dev_Colab,
    "Referencias" => "reference.md",
    "Contatos" => "contatos.md",
    # [
    #   nice_name(file) => file for
    #   file in readdir(joinpath(@__DIR__, "src")) if file != "index.md" && contains(file, ".md") && splitext(file)[2] == ".md"
    # ]...
  ],
)

deploydocs(; repo = "github.com/Mesquitaaph/MyProject.jl", push_preview = false)


# Toda vez que iniciar o REPL da Julia precisa seguir os seguintes passos:
# 1 - pkg> activate docs
# 2 - julia> include("docs/make.jl")