using Fease 
using GaussQuadrature

# ============== Funções específicas ===============

function getQuadSideNodes(Nx1::Int64, Nx2::Int64, side::Int64)::Vector{Union{Any, Int64}}
    nodes = []
    
    if side == 1
        nodes = reduce(vcat, 1:Nx1+1)
    end
    if side == 2
        nodes = reduce(vcat, Nx1+1:Nx1+1:(Nx2+1)*(Nx1+1))
    end
    if side == 3
        nodes = reduce(vcat, Nx2*(Nx1+1)+1:(Nx2+1)*(Nx1+1)) 
    end
    if side == 4
        nodes = reduce(vcat, 1:Nx1+1:Nx2*(Nx1+1)+1)
    end

    return nodes    
end

function getQuadSideElements(Nx1::Int64, Nx2::Int64, side::Int64)::Vector{Union{Any, Int64}}
    elements = []

    if side == 1
        elements = reduce(vcat, 1:Nx1)
    end
    if side == 2
        elements = reduce(vcat, Nx1:Nx1:(Nx2)*(Nx1))
    end
    if side == 3
        elements = reduce(vcat, Nx1*(Nx2-1)+1:1:(Nx2)*(Nx1)) 
    end
    if side == 4
        elements = reduce(vcat, 1:Nx1:Nx1*(Nx2-1)+1)
    end

    return elements  
end

function getSidefromNode(node::Int64, Nx1::Int64, Nx2::Int64)
    
    sides = []
    
    if node in 1:Nx1+1
        append!(sides, 1)
    end
        
    if node in Nx1+1:Nx1+1:(Nx2+1)*(Nx1+1)
        append!(sides, 2)
    end
    
    if node in Nx2*(Nx1+1)+1:(Nx2+1)*(Nx1+1)
        append!(sides, 3)
    end
    
    if node in 1:Nx1+1:Nx2*(Nx1+1)+1
        append!(sides, 4)
    end
    
    return sides
    
end

function getBoundariesOfElement(e::Int64, Nx1::Int64, Nx2::Int64, LG)
    boundaries = []

    for node in LG[:, e]
        append!(boundaries, getSidefromNode(node, Nx1, Nx2))
    end

    return unique(boundaries)
end

function getQuad(x, y, LG, elem)
    
    idx00 = LG[elem, 1]    
    idx01 = LG[elem, 2]    
    idx11 = LG[elem, 3]    
    idx10 = LG[elem, 4]
    
    X = [x[idx00], x[idx10], x[idx11], x[idx01]]
    Y = [y[idx00], y[idx10], y[idx11], y[idx01]]
    
    return X, Y
end

function drawQuadInFig(fig, xVertices, yVertices)
    for i in 1:3
        xPair = [xVertices[i], xVertices[i+1]]
        yPair = [yVertices[i], yVertices[i+1]]
        linesegments!(fig[1,1], xPair, yPair, color=:navy)
    end
    
    xPair = [xVertices[4], xVertices[1]]
    yPair = [yVertices[4], yVertices[1]]
    linesegments!(fig[1,1], xPair, yPair, color=:navy)
    
    return fig
end

function drawGrid(xValues, yValues, numX, numY, LG, windowX, windowY)
    LGᵗ = transpose(LG)
    
    f = Figure(size=(500, 500))
    Axis(f[1, 1], limits=(0, windowX, 0, windowY))
    
    for i in 1:numY
        for j in 1:numX
            
            elem = (i-1)*(numX) + j
            
            xVertices, yVertices = getQuad(xValues, yValues, LGᵗ, elem)
            drawQuadInFig(f, xVertices, yVertices)
        end
    end
    
    return f
end

function drawGridFromFile(elemX, elemY, coordsFileName, LGFileName, xWindowLim, yWindowLim)
    malhaFile = readdlm(coordsFileName)

    X, Y = malhaFile[:, 2], malhaFile[:, 3]

    LG = readdlm(LGFileName)[:, 1:4]
    LG = convert(Matrix{Int}, LG) 
    LG = LG .+ 1 # Indices 1-based (julia)
    
    return drawGrid(X, Y, elemX, elemY, transpose(LG), xWindowLim, yWindowLim)
end 

function mapper_to_x_generic(Xᵉ_a::Vector{Any}, ξ₁::Float64, ξ₂::Float64)::Vector{Float64}
    return ϕ_2D(ξ₁, ξ₂)
end

function bound_expr(num_fronteira, nos_fronteira, e_fronteira, malha, P)
    coords_nos_front = []
    
    for no in nos_fronteira
        append!(coords_nos_front, [[malha.coords[1][no], malha.coords[2][no]]])
    end

    for e in e_fronteira
        Xᵉ = elem_coords(malha, e)[2]
        
        if(num_fronteira == 1)
            for ξ in P
                x1_ξ1 = mapper_to_x_generic(Xᵉ[1], ξ, -1.0)
                x2_ξ2 = mapper_to_x_generic(Xᵉ[2], ξ, -1.0)  
                append!(coords_nos_front, [[x1_ξ1, x2_ξ2]])
            end
        end

        if(num_fronteira == 2)
            for ξ in P
                x1_ξ1 = mapper_to_x_generic(Xᵉ[1], 1.0, ξ)
                x2_ξ2 = mapper_to_x_generic(Xᵉ[2], 1.0, ξ)
                append!(coords_nos_front, [[x1_ξ1, x2_ξ2]])
            end
        end

        if(num_fronteira == 3)
            for ξ in P
                x1_ξ1 = mapper_to_x_generic(Xᵉ[1], ξ, 1.0)
                x2_ξ2 = mapper_to_x_generic(Xᵉ[2], ξ, 1.0)
                append!(coords_nos_front, [[x1_ξ1, x2_ξ2]])
            end
        end

        if(num_fronteira == 4)
            for ξ in P
                x1_ξ1 = mapper_to_x_generic(Xᵉ[1], -1.0, ξ)
                x2_ξ2 = mapper_to_x_generic(Xᵉ[2], -1.0, ξ)
                append!(coords_nos_front, [[x1_ξ1, x2_ξ2]])
            end
        end
    end

    #println(coords_nos_front)
    
    expr =  function(coords_no)
                if(coords_no in coords_nos_front)
                    return true
                else
                    return false
                end
            end
    
    return expr
end

function build_f(pontos_fronteira::Array, nos_fronteira, e_fronteira, malha, P, τ::Float64)  
    p1 = pontos_fronteira[1]
    p2 = pontos_fronteira[2]
    p3 = pontos_fronteira[3]
    p4 = pontos_fronteira[4]
    
    dy = p4[2] - p1[2]
    dx = p4[1] - p1[1]
    dl = sqrt(dx^2 + dy^2)
    _sen = dy/dl
    _cos = dx/dl

    bound1 = bound_expr(1, nos_fronteira[1], e_fronteira[1], malha, P)
    bound2 = bound_expr(2, nos_fronteira[2], e_fronteira[2], malha, P)
    bound3 = bound_expr(3, nos_fronteira[3], e_fronteira[3], malha, P)
    bound4 = bound_expr(4, nos_fronteira[4], e_fronteira[4], malha, P)
    
    f = function(ponto)
            tracao = zeros(2) 
            #println("Checking for ", ponto)
        
            if(bound1(ponto))
                tracao[1] += -τ
                tracao[2] += 0.0
            end
        
            if(bound2(ponto))
                tracao[1] += τ * _cos
                tracao[2] += τ * _sen
            end
        
            if(bound3(ponto))
                tracao[1] += τ
                tracao[2] += 0.0
            end
        
            if(bound4(ponto))
                tracao[1] += -τ * _cos
                tracao[2] += -τ * _sen
            end
        
            return tracao
        end
    
    return f 
end

# ============== Início script =====================

# Infos Elementos Finitos
Nx1 = 20
Nx2 = 20

# Infos ALI
λ₁ = 1.0
λ₂ = 1.0

n_passos = 8
Δτ = 0.5/n_passos

s1 =  1.0
s2 = -0.1
β  = 10.0^4
params = [s1, s2, β]

# Infos Display
winX = 2.0
winY = 2.0

ponto_inf_esq = (0.0, 0.0)
ponto_inf_dir = (λ₁ , 0.0)
ponto_sup_dir = (λ₁ , λ₂ )
ponto_sup_esq = (0.0, λ₂ )

nos_fronteiras = [
    getQuadSideNodes(Nx1, Nx2, 1), 
    getQuadSideNodes(Nx1, Nx2, 2), 
    getQuadSideNodes(Nx1, Nx2, 3),
    getQuadSideNodes(Nx1, Nx2, 4)
]

# Neste exemplo, do quadrilatero apoiado no chão, estas são as prescrições
presc_x1 = [1]
presc_x2 = nos_fronteiras[1]
presc = [presc_x1, presc_x2]

baseType = BaseTypes.linearLagrange
n_dir = 2
malha = monta_malha_2D_uniforme(baseType, Nx1, Nx2, ponto_inf_esq, ponto_sup_dir, n_dir, presc)

pts_fronteira = [ponto_inf_esq, ponto_inf_dir, ponto_sup_dir, ponto_sup_esq]

elementos_fronteiras = [
    getQuadSideElements(Nx1, Nx2, 1),
    getQuadSideElements(Nx1, Nx2, 2),
    getQuadSideElements(Nx1, Nx2, 3),
    getQuadSideElements(Nx1, Nx2, 4),
]

P, W = legendre(2)
τ = 1 * Δτ
f = build_f(pts_fronteira, nos_fronteiras, elementos_fronteiras, malha, P, τ) 

# coords = joinToPoints(X₁, X₂)
# nodes = (unique(Iterators.flatten(presc)))
# values = repeat([[0.0; 0.0]], size(nodes)[1])

# g = build_g(nodes, coords, values)    

# s1, s2, β =  params
# h  = -s2/s1
# # Atualizações de p?? Tr(H) ≈ 0 sempre, pode ficar de fora
# p = 1-h
# params = [s1, s2, p, β]

# ne = Nx1 * Nx2
# npts = (Nx1 + 1) * (Nx2 + 1)

# # F e T iniciais: Identidade e Matriz nula
# F_def = repeat([Matrix(1.0I, 2, 2)], ne)
# T = repeat([Matrix(0.0I, 2, 2)], ne)

# #X_novo = []
# #K_vec = []
# #F_vec = []
# sol_vec = []
# infos_primarias = []
# infos_secundarias = []
# for iter in 1:n_passos
    
#     K, F = monta_K_F_global(params, f, pts_fronteira, g, presc, Nx1, Nx2, X, F_def, T, m, EQ, LG)
#     c = K\F
    
#     # u em t+1
#     X_novo = soma_solução(X, c, npts, m, EQ)
#     #display(drawGrid(X_novo[:, 1], X_novo[:, 2], Nx1, Nx2, LG, winX, winY))
    
#     novo_F, novo_T = atualiza_ALI(X, X_novo, F_def, T, ne, LG)
    
#     F_def = copy(novo_F)
#     T = copy(novo_T)
#     X = copy(X_novo)
#     append!(sol_vec, [X])

#     infos1, infos2 = calcula_erro(iter, F_def, s1, s2, τ, pts_fronteira, ne)
#     append!(infos_primarias, [infos1])
#     append!(infos_secundarias, [infos2])

#     # Fim do passo atual, prepara pra novo passo
#     pts_fronteira = [X[first(fronteiras[1]), :], 
#                         X[last(fronteiras[1]), :], 
#                         X[last(fronteiras[3]), :], 
#                         X[first(fronteiras[3]), :]]::Vector{Vector{Float64}}
    
#     τ = Δτ * (iter + 1)
#     f = build_f(pts_fronteira, fronteiras, e_fronteiras, X, P, LG, τ)  #function(x) return [0.0, 0.0] end
# end

# return sol_vec, infos_primarias, infos_secundarias

# function calcL(s1::Float64, s2::Float64, β::Float64, B, inv_B)
#     L = zeros(2,2,2,2) 
#     delta = 1.0I
    
#     for i in 1:2
#         for j in 1:2
#             for k in 1:2
#                 for l in 1:2
#                     L[i, j, k, l] = β*delta[k, l]*delta[i, j] + s1*(delta[i, k]*B[l, j] + B[i, l]*delta[j, k]) -
#                                     s2*(inv_B[i, k] * delta[l, j] + delta[l, i]*inv_B[k, j])
#                 end
#             end
#         end
#     end
    
#     return L
# end

# function build_g(nodes, coords, values)::Function
#     nodes_coords = coords[nodes, :]
    
#     g = function(x) 
#             if(x in nodes_coords)
#                 return values[findfirst(item -> item == x, nodes_coords)]
#             else
#                 return [0.0; 0.0]
#             end
#         end
    
#     return g    
# end

# function build_f(pontos_fronteira::Array, nos_fronteira, e_fronteira, X, P, LG, τ::Float64)  
#     p1 = pontos_fronteira[1]
#     p2 = pontos_fronteira[2]
#     p3 = pontos_fronteira[3]
#     p4 = pontos_fronteira[4]
    
#     dy = p4[2] - p1[2]
#     dx = p4[1] - p1[1]
#     dl = sqrt(dx^2 + dy^2)
#     _sen = dy/dl
#     _cos = dx/dl

#     bound1 = bound_expr(1, nos_fronteira[1], e_fronteira[1], X, P, LG)
#     bound2 = bound_expr(2, nos_fronteira[2], e_fronteira[2], X, P, LG)
#     bound3 = bound_expr(3, nos_fronteira[3], e_fronteira[3], X, P, LG)
#     bound4 = bound_expr(4, nos_fronteira[4], e_fronteira[4], X, P, LG)
    
#     f = function(ponto)
#             tracao = zeros(2) 
#             #println("Checking for ", ponto)
        
#             if(bound1(ponto))
#                 tracao[1] += -τ
#                 tracao[2] += 0.0
#             end
        
#             if(bound2(ponto))
#                 tracao[1] += τ * _cos
#                 tracao[2] += τ * _sen
#             end
        
#             if(bound3(ponto))
#                 tracao[1] += τ
#                 tracao[2] += 0.0
#             end
        
#             if(bound4(ponto))
#                 tracao[1] += -τ * _cos
#                 tracao[2] += -τ * _sen
#             end
        
#             return tracao
#         end
    
#     return f 
# end

