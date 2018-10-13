__precompile__(true)

"""
Module `Approximations.jl` -- polygonal approximation of convex sets through
support vectors.
"""
module Approximations

using LazySets, Requires

export approximate,
       ballinf_approximation,
       box_approximation,
       decompose,
       diameter,
       norm,
       overapproximate,
       radius,
       symmetric_interval_hull,
       BoxDirections,
       BoxDiagDirections,
       OctDirections

include("../compat.jl")

const TOL(N::Type{Float64}) = eps(N)
const TOL(N::Type{Float32}) = eps(N)
const TOL(N::Type{Rational{INNER}}) where {INNER} = zero(N)

const DIR_EAST(N) = [one(N), zero(N)]
const DIR_NORTH(N) = [zero(N), one(N)]
const DIR_WEST(N) = [-one(N), zero(N)]
const DIR_SOUTH(N) = [zero(N), -one(N)]

include("iterative_refinement.jl")
include("box_approximations.jl")
include("template_directions.jl")
include("overapproximate.jl")
include("decompositions.jl")
#include("upper_bounds.jl")

end # module
