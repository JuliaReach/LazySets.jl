"""
Module `Approximations.jl` -- polygonal approximation of sets.
"""
module Approximations

import Base: convert
import ..LazySets: dim, project, □

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
       isbounding,
       overapproximate_norm

import IntervalArithmetic as IA

using LinearAlgebra: /, I, dot, norm, normalize, nullspace, transpose
using SparseArrays: SparseVector, sparsevec, spzeros
using Requires: @require

using ReachabilityBase.Arrays: At_mul_B, SingleEntryVector, rectify,
                               remove_zero_columns, uniform_partition
using ReachabilityBase.Comparison: _isapprox, _leq, _geq, _rtol, isapproxzero
using ReachabilityBase.Subtypes: subtypes

using ..LazySets
using ..LazySets: ACS, default_lp_solver, _isbounded_stiemke, require, linprog,
                  is_lp_optimal, _normal_Vector, default_sdp_solver,
                  get_exponential_backend, _expmv, second, @assert, _box_radius
using ..LazySets.JuMP: Model, set_silent, @variable, @constraint, optimize!,
                       value, @NLobjective, @objective
using ..LazySets.MatrixZonotopeModule: _rowwise_zonotope_norm

include("box_approximation.jl")
include("iterative_refinement.jl")
include("symmetric_interval_hull.jl")
include("ballinf_approximation.jl")
include("template_directions.jl")
include("overapproximate.jl")
include("overapproximate_interval.jl")
include("overapproximate_zonotope.jl")
include("overapproximate_cartesianproductarray.jl")
include("underapproximate.jl")
include("approximate.jl")
include("decompositions.jl")
include("hausdorff_distance.jl")
include("overapproximate_norm.jl")
include("overapproximate_expmap.jl")
include("overapproximate_matrixzonotope.jl")
include("init.jl")

end # module
