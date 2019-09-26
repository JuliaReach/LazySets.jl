__precompile__(true)

"""
Module `Parallel.jl` -- LazySets algorithms that are parallelized.
"""
module Parallel

using LazySets
using ..Assertions: @assert, activate_assertions
# activate assertions by default
activate_assertions(Parallel)

using SharedArrays: SharedMatrix, SharedVector, indexpids
using Distributed: remotecall_wait, procs

#=======================================================
Utility functions for distribution of tasks in parallel
=======================================================#
include("distribute.jl")

#==================================================
Approximations using boxes implemented in parallel
==================================================#
include("box_approximations.jl")

end # module
