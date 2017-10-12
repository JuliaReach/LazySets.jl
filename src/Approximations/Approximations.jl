__precompile__(true)

"""
Module `Approximations.jl` -- polygonal approximation of convex sets through
support vectors.
"""
module Approximations

using LazySets

export decompose, overapproximate, approximate,
       radius_approximation, diameter_approximation, box_approximation,
       ballinf_approximation, symmetric_interval_hull

const TOL_DIR = 1e-6
const DIR_EAST = [1., 0.]
const DIR_NORTH = [0., 1.]
const DIR_WEST = [-1., 0.]
const DIR_SOUTH = [0., -1.]

include("iterative_refinement.jl")
include("box_approximations.jl")
include("decompositions.jl")

end # module
