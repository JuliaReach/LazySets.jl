__precompile__(true)

"""
Module `Arrays.jl` -- Auxiliary machinery for vectors and matrices.
"""
module Arrays

using LinearAlgebra, SparseArrays
using ..Assertions: @assert, activate_assertions
# activate assertions by default
activate_assertions(Arrays)

include("matrix_operations.jl")
include("vector_operations.jl")
include("matrix_vector_operations.jl")
include("SingleEntryVector.jl")

end # module
