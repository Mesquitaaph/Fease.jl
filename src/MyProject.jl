module MyProject

export format_num, measure_func

export Malha, monta_malha

export examples, RunValues

export BaseTypes, LocalBases, monta_base

export solveSys, PHI, PHI_original

export convergence_test!

export single_run


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

include("test.jl")

end
