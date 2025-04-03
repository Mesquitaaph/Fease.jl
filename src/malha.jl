struct Malha
  base
  ne::Int64
  neq::Int64
  X
  dx
  EQ
  LG
  EQoLG::Matrix{Int64}
  a
  b
end

function montaLG_geral(Nx1::Int64, Nx2::Int64 = 0)
  n_dim = 1 + (Nx2 == 0 ? 0 : 1)
  ne = Nx1 * (Nx2 == 0 ? 1 : Nx2)

  LG = zeros(Int64, 2^n_dim, ne)

  # Monta o primeiro bloco
  LG1 = [i for i in 1:Nx1]
  for e in 1:Nx1
    for a in 1:2^n_dim
      @inbounds LG[a,e] = LG1[e] + ((a-1) % 2) + (Nx1 + 1)*floor((a-1)/2)
    end
  end

  # Preenche o restante dos blocos baseado no primeiro
  for e in Nx1+1:ne
    @inbounds LG[:,e] .= LG[:,e-Nx1] .+ (Nx1+1)
  end

  return LG
end

function montaLG(ne, base)
  LG = zeros(Int64, ne, base.nB)

  LG1 = [(base.nB-1)*i - (base.nB-2) for i in 1:ne]
  
  for e in 1:ne
      for a in 1:base.nB
          @inbounds LG[e,a] = LG1[e] + (a-1)
      end
  end

  return LG'
end

function monta_LG_2D(Nx1::Int64, Nx2::Int64)::Matrix{Int64}
  # Define o número de funções φ no eixo x₁ e x₂
  nx1 = Nx1 + 1
  nx2 = Nx2 + 1

  # M[:,j] contém a numeração da primeira linha do "Bloco j"
  M = (1:nx1-1) .+ (0:nx1:(nx2-2)*nx1)'

  # LG[1,:] contém a numeração global da primeira função local de cada elemento 
  linha1 = reshape(M, 1, :)

  # Constrói a matriz LG
  LG = vcat(linha1, linha1 .+ 1,  linha1 .+ nx1, linha1 .+ (nx1+1))

  return LG
end

function montaEQ_geral(Nx1::Int64, Nx2::Int64 = 0)
  n_dim = 1 + (Nx2 == 0 ? 0 : 1)
  ne = Nx1 * (Nx2 == 0 ? 1 : Nx2)

  neq = (Nx1-1) * (Nx2 == 0 ? 1 : Nx2-1)

  EQ = fill(neq+1, (Nx1+1) * (Nx2 == 0 ? 1 : Nx2+1))
  
  contador = 1
  for y in 1:Nx2+1
    for x in 1:Nx1+1
      if !(x == 1 || x == Nx1+1 || (Nx2 > 0 && (y == 1 || y == Nx2+1)))
        idx = x + (y-1)*(Nx1+1)
        EQ[idx] = contador
        contador += 1
      end
    end
  end

  return EQ
end

function montaEQ(ne, neq, base)
  EQ = zeros(Int64, (base.nB-1)*ne+1)

  EQ[1] = neq+1; EQ[end] = neq+1
  
  for i in 2:(base.nB-1)*ne
      @inbounds EQ[i] = i-1
  end

  return EQ
end

function monta_EQ_2D(Nx1::Int64, Nx2::Int64)
  # Define o número de funções φ no eixo x₁ e x₂
  nx1 = Nx1 + 1
  nx2 = Nx2 + 1
  
  # Calcula o número de funções globais φ que compõem a base do espaço Vₘ
  m = (nx1-2) * (nx2-2)

  # Inicializa o vetor EQ preenchido com m+1
  EQ = fill(m+1, nx1 * nx2)

  # Vetor contendo os índices das funções globais φ que compõem a base do espaço Vₘ
  L = reshape( (0:nx1-3) .+ (nx1+2:nx1:(nx2-2)*nx1+2)' , :,1)

  # Atribui os valores de 1 até m as funções globais φ que compõem a base do espaço Vₘ
  EQ[L] = 1:m

  return EQ
end

function monta_malha(ne, base, a, b)
  dx = 1/ne;
  X = a:dx:b

  neq = base.neq

  EQ = montaEQ(ne, neq, base); LG = montaLG(ne, base)
  EQoLG = EQ[LG]

  return Malha(base, ne, neq, X, dx, EQ, LG, EQoLG, a, b)
end