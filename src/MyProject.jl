module MyProject

export format_num, measure_func

export examples, RunValues

export BaseTypes, LocalBases

export solveSys, PHI, PHI_original

export convergence_test!

export single_run


# Dependencias do modulo
using 
  BenchmarkTools, GaussQuadrature, Plots, BandedMatrices, Printf, 
  DataFrames, Latexify, Statistics, SparseArrays


include("utils.jl")

include("examples.jl")

include("bases.jl")

include("serial.jl")
# include("FEM2D.jl")

include("convergence.jl")

include("test.jl")

end
