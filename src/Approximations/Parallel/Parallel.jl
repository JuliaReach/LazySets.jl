__precompile__(true)

"""
Module `Approximations.Parallel.jl` -- parallel implementations of some approximation
algorithms.
"""
module Parallel

using LazySets

#=======================================================
Utility functions for distribution of tasks in parallel
=======================================================#
include("distribute.jl")

#==================================================
Approximations using boxes implemented in parallel
==================================================#
include("box_approximation.jl")

end # module