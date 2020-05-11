__precompile__(true)

"""
Module `Approximations.jl` -- polygonal approximation of convex sets through
support vectors.
"""
module Approximations

using LazySets, LazySets.Arrays, Requires, LinearAlgebra, SparseArrays,
      MathProgBase
using LazySets: _isapprox, _rtol, _normal_Vector, isapproxzero, default_lp_solver
using ..Assertions: @assert, activate_assertions
# activate assertions by default
activate_assertions(Approximations)

export approximate,
       project,
       ballinf_approximation,
       box_approximation, interval_hull,
       decompose,
       overapproximate,
       underapproximate,
       box_approximation_symmetric, symmetric_interval_hull,
       BoxDirections,
       BoxDiagDirections,
       OctDirections,
       PolarDirections,
       SphericalDirections,
       CustomDirections,
       isbounding

const DIR_EAST(N) = [one(N), zero(N)]
const DIR_NORTH(N) = [zero(N), one(N)]
const DIR_WEST(N) = [-one(N), zero(N)]
const DIR_SOUTH(N) = [zero(N), -one(N)]

include("iterative_refinement.jl")
include("box_approximations.jl")
include("template_directions.jl")
include("overapproximate.jl")
include("underapproximate.jl")
include("decompositions.jl")
include("hausdorff_distance.jl")
include("init.jl")

end # module
