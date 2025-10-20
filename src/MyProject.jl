module MyProject

# Dependencias do modulo
using
      BenchmarkTools, GaussQuadrature, Plots, BandedMatrices, Printf,
      DataFrames, Latexify, Statistics, SparseArrays, LinearAlgebra, Random,
      ProfileView

include("utils.jl")

include("bases.jl")

include("malha.jl")

include("examples.jl")

include("phis.jl")

include("serial.jl")
include("FEM2D.jl")
include("generalizado.jl")

include("convergence.jl")

include("plots.jl")

include("testes/include_testes.jl")

export format_num, measure_func, test_revise

export Malha, monta_malha_1D_uniforme, monta_malha_2D_uniforme, malha2D_adiciona_ruido

export examples_1D, examples_2D, RunValues, TermosEquacao

export BaseTypes, LocalBases, monta_base

export montaK_geral, montaF_geral, solveSys_geral, solve_sys, solve_sys_poisson,
       monta_u_aproximada

export dot

export plot_malha_2D, plot_solucao_aproximada

export single_run_1D, single_run_2D

export erro_L2, convergence_test!

using JuliaFormatter
export format

end
