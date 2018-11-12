__precompile__(true)

"""
Module `Parallel.jl` -- LazySets algorithms that are parallelized.
"""
module Parallel

using LazySets

# compatibility of different julia versions
@static if VERSION >= v"0.7-"
    using SharedArrays: SharedMatrix, SharedVector, indexpids
    using Distributed: remotecall_wait, procs
end

#=======================================================
Utility functions for distribution of tasks in parallel
=======================================================#
include("distribute.jl")

#==================================================
Approximations using boxes implemented in parallel
==================================================#
include("box_approximations.jl")

end # module
