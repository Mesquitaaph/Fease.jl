using Documenter, Fease

DocMeta.setdocmeta!(Fease, :DocTestSetup, :(using Fease); recursive = true)

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

Tutoriais_Desenvolvimento = "De Desenvolvimento" => [
  "Novo termo de Equação" => "tutoriais/novo-termo.md",
]

Tutoriais = "Tutoriais" => [
  "Tutoriais" => "tutoriais/index.md",
  Tutoriais_Desenvolvimento
]

makedocs(;
  modules = [Fease],
  doctest = true,
  linkcheck = false, # Rely on Lint.yml/lychee for the links
  authors = "Raphael Mesquita <raphaelfcm@ic.ufrj.br> and contributors",
  repo = "https://github.com/Mesquitaaph/Fease.jl/blob/{commit}{path}#{line}",
  sitename = "Fease.jl",
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical = "https://Mesquitaaph.github.io/Fease.jl",
    # assets = ["assets/style.css"],
  ),
  pages = [
    "Início" => "index.md",
    Para_Comecar,
    "Método de Elementos Finitos" => "metodo-elementos-finitos.md",
    Tutoriais,
    Dev_Colab,
    "Referencias" => "reference.md",
    "Contatos" => "contatos.md",
    # [
    #   nice_name(file) => file for
    #   file in readdir(joinpath(@__DIR__, "src")) if file != "index.md" && contains(file, ".md") && splitext(file)[2] == ".md"
    # ]...
  ],
)

deploydocs(; repo = "github.com/Mesquitaaph/Fease.jl", push_preview = false)


# Toda vez que iniciar o REPL da Julia precisa seguir os seguintes passos:
# 1 - pkg> activate docs
# 2 - julia> include("docs/make.jl")