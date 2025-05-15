module MyProject

# Dependencias do modulo
using 
  BenchmarkTools, GaussQuadrature, Plots, BandedMatrices, Printf, 
  DataFrames, Latexify, Statistics, SparseArrays, LinearAlgebra, Random

include("utils.jl")

include("bases.jl")

include("malha.jl")

include("examples.jl")

include("phis.jl")

include("serial.jl")
include("FEM2D.jl")
include("generalizado.jl")

include("convergence.jl")

include("testes/testes_malha.jl")
include("testes/testes_simulacao.jl")


export format_num, measure_func

export Malha, monta_malha_1D_uniforme, monta_malha_2D_uniforme, malha2D_adiciona_ruido

export examples, RunValues

export BaseTypes, LocalBases, monta_base

export solveSys, PHI, PHI_original

export exemplo1

export convergence_test!

export solveSys_geral, montaK_geral, montaF_geral

export dot

end
