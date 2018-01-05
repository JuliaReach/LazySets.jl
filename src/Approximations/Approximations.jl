__precompile__(true)

"""
Module `Approximations.jl` -- polygonal approximation of convex sets through
support vectors.
"""
module Approximations

using LazySets

export approximate,
       ballinf_approximation,
       box_approximation,
       decompose,
       diameter,
       norm,
       overapproximate,
       radius,
       symmetric_interval_hull

const TOL = eps(Float64) # TODO: use eps(N)

const DIR_EAST = [1., 0.] # TODO: use given precision N
const DIR_NORTH = [0., 1.]
const DIR_WEST = [-1., 0.]
const DIR_SOUTH = [0., -1.]

include("iterative_refinement.jl")
include("box_approximations.jl")
include("decompositions.jl")

end # module
