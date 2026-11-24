using Fease
using Plots

function caso_2D()
    # Define o número de sub-intervalos no eixo horizontal e vertical, respectivamente.
    Nx1, Nx2 = 4, 3

    # Define o tipo da base de funções interpoladoras do subespaço aproximado Vₘ.
    baseType = BaseTypes.linearLagrange

    # Define extremos do intervalo da malha.
    a = (0.0, 0.0)
    b = (1.0, 1.0)

    # Define nos prescritos (neste caso, todo o contorno)
    nos_prescritos = [
        getQuadSideNodes(Nx1, Nx2, 1),
        getQuadSideNodes(Nx1, Nx2, 2),
        getQuadSideNodes(Nx1, Nx2, 3),
        getQuadSideNodes(Nx1, Nx2, 4)
    ]

    # Formata para Vector{Vector{Any}}, para ficar compatível com os algoritmos
    nos_prescritos = [collect(Iterators.flatten(nos_prescritos))]

    # Monta a fronteira
    fronteira = Fease.monta_fronteira_2D_uniforme(a, b, Nx1, Nx2, nos_prescritos)

    # Graus de liberdade do problema(de quantas variáveis independentes depende a solução u)
    n_graus_liberdade = 1

    # Define a malha com os valores atribuídos acima.
    malha = Fease.monta_malha_2D_uniforme(baseType, Nx1, Nx2, n_graus_liberdade, fronteira)

    # Define os parâmetros da equação a ser resolvida.
    α = 1.0
    β = 1.0
    f = (x₁, x₂) -> (2 * α * π^2 + β) * sin(π * x₁) * sin(π * x₂)

    # Define o pseudo operador linear a(u,v).
    function pseudo_a(termos_equacao)
        (; ∇u, ∇v, u, v) = termos_equacao

        return β * dot(u, v) + α * dot(∇u, ∇v)
    end

    # Monta e resolve o sistema linear relacionado a esta equação.
    c = solve_sys(f, malha, pseudo_a)

    # Imprime solução aproximada.
    show(c)

    # Plota solução aproximada.
    plot_solucao_aproximada(c, malha, false)

    println("\n Feche o plot e pressione Enter pra continuar...")
    readline()

    # Estudo da convergencia do erro

    # Define função que constrói uma malha que depende só do número de sub-intervalos nos dois eixos.
    function monta_malha(NX)
        # Define nos prescritos (neste caso, todo o contorno)
        nos_prescritos = [
            getQuadSideNodes(NX[1], NX[2], 1),
            getQuadSideNodes(NX[1], NX[2], 2),
            getQuadSideNodes(NX[1], NX[2], 3),
            getQuadSideNodes(NX[1], NX[2], 4)
        ]

        # Formata para Vector{Vector{Any}}, para ficar compatível com os algoritmos
        nos_prescritos = [collect(Iterators.flatten(nos_prescritos))]

        # Monta a fronteira
        fronteira = Fease.monta_fronteira_2D_uniforme(a, b, NX[1], NX[2], nos_prescritos)

        # Graus de liberdade do problema(de quantas variáveis independentes depende a solução u)
        n_graus_liberdade = 1

        # Define a malha com os valores atribuídos acima.
        return Fease.monta_malha_2D_uniforme(baseType, NX[1], NX[2], n_graus_liberdade, fronteira)
    end

    # Estudo será realizado com número de elementos em potências de 2.
    errsize = 7
    NE = 2 .^ [2:1:errsize;]
    H = 1 ./ NE
    E = zeros(length(NE))

    # Define a solução analítica da equação.
    u = (x₁, x₂) -> sin(π * x₁) * sin(π * x₂)

    # Calcula o erro para cada uma das quantidades de elementos finitos.
    n_dim = 2
    convergence_test!(E, NE, n_dim, monta_malha, pseudo_a, f, u)

    # Plota o resultado do estudo
    Plots.plot(H, E, xaxis=:log10, yaxis=:log10, label="Erro")
    return Plots.plot!(H, H .^ monta_base(baseType, 2).nB, xaxis=:log10,
        yaxis=:log10, label="H²")
end

# Chama a função definida acima.
caso_2D()