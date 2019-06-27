__precompile__(true)

"""
Module `Arrays.jl` -- Auxiliary machinery for vectors and matrices.
"""
module Arrays

using LazySets, LinearAlgebra, SparseArrays

include("matrix_operations.jl")
include("vector_operations.jl")
include("matrix_vector_operations.jl")
include("SingleEntryVector.jl")

end # module
