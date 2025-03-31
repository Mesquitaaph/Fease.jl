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

function montaEQ(ne, neq, base)
  EQ = zeros(Int64, (base.nB-1)*ne+1, 1)

  EQ[1] = neq+1; EQ[end] = neq+1
  
  for i in 2:(base.nB-1)*ne
      @inbounds EQ[i] = i-1
  end

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