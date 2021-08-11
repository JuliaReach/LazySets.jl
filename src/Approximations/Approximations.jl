__precompile__(true)

"""
Module `Approximations.jl` -- polygonal approximation of convex sets through
support vectors.
"""
module Approximations

using LazySets, LazySets.Arrays, Requires, LinearAlgebra, SparseArrays,
      MathProgBase

using LazySets: _isapprox, _leq, _geq, _rtol, _normal_Vector, isapproxzero,
                default_lp_solver, _isbounded_stiemke, require

import LazySets: project

using ..Assertions: @assert, activate_assertions
# activate assertions by default
activate_assertions(Approximations)

export approximate,
       ballinf_approximation,
       box_approximation, interval_hull,
       decompose,
       overapproximate,
       underapproximate,
       box_approximation_symmetric, symmetric_interval_hull,
       BoxDirections,
       DiagDirections,
       BoxDiagDirections,
       OctDirections,
       PolarDirections,
       SphericalDirections,
       CustomDirections,
       isbounding

include("box_approximation.jl")
include("iterative_refinement.jl")
include("symmetric_interval_hull.jl")
include("ballinf_approximation.jl")
include("template_directions.jl")
include("overapproximate.jl")
include("underapproximate.jl")
include("approximate.jl")
include("decompositions.jl")
include("hausdorff_distance.jl")
include("init.jl")

end # module
