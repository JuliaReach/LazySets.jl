using BenchmarkTools, LazySets
using LinearAlgebra, SparseArrays, Random

SUITE = BenchmarkGroup()  # parent BenchmarkGroup to contain our suite

include("supports.jl")
include("linear_map.jl")
