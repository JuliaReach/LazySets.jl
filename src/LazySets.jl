module LazySets

using Reexport

include("API/API.jl")
@reexport using .API
import .API: eltype, extrema, isdisjoint, isempty, \, ∈, ≈, ==, ⊆,
             rand, norm, permute, distance, rectify,
             affine_map, an_element, area, center, complement, concretize,
             constraints_list, constraints, convex_hull, diameter, dim, exponential_map,
             high, is_interior_point, ispolyhedral, isbounded, isboundedtype,
             isconvextype, isempty, isoperation, isoperationtype, isuniversal,
             linear_map, low, norm, project, radius, reflect, sample, scale,
             scale!, support_function, ρ, support_vector, σ, surface, translate,
             translate!, vertices_list, vertices, volume,
             cartesian_product, difference, distance, exact_sum, intersection,
             is_intersection_empty, isequivalent, ⊂, linear_combination,
             minkowski_difference, pontryagin_difference, minkowski_sum

@reexport import LinearAlgebra: ×, normalize, normalize!
import Base: copy, rationalize, \
import RecipesBase: apply_recipe
export Arrays, subtypes

using LinearAlgebra, RecipesBase, Requires, SparseArrays
import GLPK, JuMP, Random, ReachabilityBase
import IntervalArithmetic as IA
using LinearAlgebra: checksquare
using Random: AbstractRNG, GLOBAL_RNG, SamplerType, randperm

@static if VERSION < v"1.9"
    stack(vecs) = hcat(vecs...)
end

# ================
# ReachabilityBase
# ================

import ReachabilityBase.Assertions
using ReachabilityBase.Assertions: @assert
include("Utils/assertions.jl")

using ReachabilityBase.Arrays
using ReachabilityBase.Basetype
using ReachabilityBase.Commutative
using ReachabilityBase.Comparison
using ReachabilityBase.Distribution
using ReachabilityBase.Iteration
using ReachabilityBase.Require
using ReachabilityBase.Subtypes

# =================
# External packages
# =================
include("Initialization/init_GLPK.jl")
include("Initialization/init_IntervalArithmetic.jl")
include("Initialization/init_JuMP.jl")

# ===================
# Auxiliary functions
# ===================
include("Utils/numbers.jl")
include("Utils/helper_functions.jl")
include("Utils/macros.jl")
include("Utils/matrix_exponential.jl")
include("Utils/lp_solvers.jl")
include("Utils/sdp_solvers.jl")
include("Utils/file_formats.jl")

# ==================
# Abstract set types
# ==================
include("Interfaces/LazySet.jl")
include("Interfaces/ConvexSet.jl")
# include("Interfaces/AbstractStar.jl")
include("Interfaces/AbstractPolynomialZonotope.jl")
include("Interfaces/AbstractSparsePolynomialZonotope.jl")
include("Interfaces/AbstractPolyhedron.jl")

include("Sets/HalfSpace/HalfSpaceModule.jl")  # must come before AbstractPolyhedron_functions
@reexport using ..HalfSpaceModule: HalfSpace, LinearConstraint,
                                   halfspace_left, halfspace_right,
                                   iscomplement
using ..HalfSpaceModule: is_tighter_same_dir_2D, _non_element_halfspace,
                         _normalize_halfspace, _normal_Vector

include("Interfaces/AbstractPolyhedron_functions.jl")
include("Interfaces/AbstractPolytope.jl")
include("Interfaces/AbstractCentrallySymmetric.jl")
include("Interfaces/AbstractCentrallySymmetricPolytope.jl")
include("Interfaces/AbstractCentrallySymmetric_functions.jl")
include("Interfaces/AbstractZonotope.jl")
include("Interfaces/AbstractHyperrectangle.jl")
include("Interfaces/AbstractPolygon.jl")
include("Interfaces/AbstractSingleton.jl")
include("Interfaces/AbstractHPolygon.jl")
include("Interfaces/AbstractAffineMap.jl")
include("Interfaces/AbstractBallp.jl")

# =============================
# Types representing basic sets
# =============================

include("Sets/Ball1/Ball1Module.jl")
@reexport using ..Ball1Module: Ball1

include("Sets/Ball2/Ball2Module.jl")
@reexport using ..Ball2Module: Ball2

include("Sets/BallInf/BallInfModule.jl")
@reexport using ..BallInfModule: BallInf

include("Sets/Ballp/BallpModule.jl")
@reexport using ..BallpModule: Ballp

include("Sets/DensePolynomialZonotope/DensePolynomialZonotopeModule.jl")
@reexport using ..DensePolynomialZonotopeModule: DensePolynomialZonotope

include("Sets/Ellipsoid/EllipsoidModule.jl")
@reexport using ..EllipsoidModule: Ellipsoid, shape_matrix

include("Sets/EmptySet/EmptySetModule.jl")
@reexport using ..EmptySetModule: EmptySet, ∅
using ..EmptySetModule: _convex_hull_emptyset, _difference_emptyset, _difference_emptyset2,
                        _distance_emptyset, _isdisjoint_emptyset, _issubset_emptyset,
                        _issubset_emptyset2, _linear_combination_emptyset,
                        _minkowski_difference_emptyset, _minkowski_difference_emptyset2,
                        _minkowski_sum_emptyset

include("Sets/HParallelotope/HParallelotopeModule.jl")
@reexport using ..HParallelotopeModule: HParallelotope,
                                        directions,
                                        base_vertex,
                                        extremal_vertices,
                                        offset

include("Sets/HPolygon/HPolygonModule.jl")
@reexport using ..HPolygonModule: HPolygon

include("Sets/HPolygonOpt.jl")

include("Sets/HPolytope/HPolytopeModule.jl")
@reexport using ..HPolytopeModule: HPolytope

include("Sets/HPolyhedron/HPolyhedronModule.jl")  # must come after HPolytopeModule
@reexport using ..HPolyhedronModule: HPolyhedron
using ..HPolyhedronModule: HPoly

include("Sets/Hyperplane/HyperplaneModule.jl")
@reexport using ..HyperplaneModule: Hyperplane
using ..HyperplaneModule: _an_element_helper_hyperplane,
                          _constraints_list_hyperplane,
                          _isdisjoint_hyperplane_hyperplane,
                          _reflect_point_hyperplane,
                          _σ_hyperplane_halfspace

include("Sets/Hyperrectangle/HyperrectangleModule.jl")
@reexport using ..HyperrectangleModule: Hyperrectangle

include("Sets/Line/LineModule.jl")
@reexport using ..LineModule: Line, direction

include("Sets/Line2D/Line2DModule.jl")
@reexport using ..Line2DModule: Line2D
using ..Line2DModule: _linear_map_hrep_helper

include("Sets/LineSegment/LineSegmentModule.jl")
@reexport using ..LineSegmentModule: LineSegment

include("Sets/Polygon/PolygonModule.jl")
@reexport using ..PolygonModule: Polygon

include("Sets/SimpleSparsePolynomialZonotope/SimpleSparsePolynomialZonotopeModule.jl")
@reexport using ..SimpleSparsePolynomialZonotopeModule: SimpleSparsePolynomialZonotope,
                                                        SSPZ,
                                                        quadratic_map

"""
    PolynomialZonotope = SimpleSparsePolynomialZonotope

Alias for `SimpleSparsePolynomialZonotope`.
"""
const PolynomialZonotope = SimpleSparsePolynomialZonotope

include("Sets/Singleton/SingletonModule.jl")
@reexport using ..SingletonModule: Singleton

include("Sets/SparsePolynomialZonotope/SparsePolynomialZonotopeModule.jl")
@reexport using ..SparsePolynomialZonotopeModule: SparsePolynomialZonotope, SPZ,
                                                  indexvector
using ..SparsePolynomialZonotopeModule: uniqueID

include("Sets/Star/StarModule.jl")
@reexport using ..StarModule: Star, basis, predicate

include("Sets/Tetrahedron/TetrahedronModule.jl")
@reexport using ..TetrahedronModule: Tetrahedron

include("Sets/Universe/UniverseModule.jl")
@reexport using ..UniverseModule: Universe
using ..UniverseModule: _convex_hull_universe, _difference_universe, _difference_universe2,
                        _distance_universe, _intersection_universe, _isdisjoint_universe,
                        _issubset_universe, _issubset_universe2, _minkowski_sum_universe

include("Sets/VPolygon/VPolygonModule.jl")
@reexport using ..VPolygonModule: VPolygon
using ..VPolygonModule: _σ_helper

include("Sets/VPolytope/VPolytopeModule.jl")
@reexport using ..VPolytopeModule: VPolytope

include("Sets/ZeroSet/ZeroSetModule.jl")
@reexport using ..ZeroSetModule: ZeroSet

include("Sets/Zonotope/ZonotopeModule.jl")
@reexport using ..ZonotopeModule: Zonotope,
                                  remove_zero_generators,
                                  linear_map!,
                                  split!
using ..ZonotopeModule: _split

include("LazyOperations/UnionSet.jl")  # must come before IntervalModule

include("Sets/Interval/IntervalModule.jl")
@reexport using ..IntervalModule: Interval

# =================================
# Types representing set operations
# =================================
include("LazyOperations/Bloating.jl")
include("LazyOperations/CartesianProduct.jl")
include("LazyOperations/CartesianProductArray.jl")
include("LazyOperations/Complement.jl")
include("LazyOperations/ConvexHull.jl")
include("LazyOperations/ConvexHullArray.jl")
include("LazyOperations/ExponentialMap.jl")
include("LazyOperations/ExponentialProjectionMap.jl")
include("LazyOperations/Intersection.jl")
include("LazyOperations/IntersectionArray.jl")
include("LazyOperations/LinearMap.jl")
include("LazyOperations/InverseLinearMap.jl")
include("LazyOperations/AffineMap.jl")  # must come after LinearMap

# must come after AffineMap
const STAR{N,VN<:AbstractVector{N},MN<:AbstractMatrix{N},PT<:AbstractPolyhedron{N}} = AffineMap{N,
                                                                                                PT,
                                                                                                N,
                                                                                                MN,
                                                                                                VN}

include("LazyOperations/MinkowskiSum.jl")
include("LazyOperations/MinkowskiSumArray.jl")
include("LazyOperations/CachedMinkowskiSumArray.jl")
include("LazyOperations/QuadraticMap.jl")
include("LazyOperations/ResetMap.jl")
include("LazyOperations/SymmetricIntervalHull.jl")
include("LazyOperations/Translation.jl")
include("LazyOperations/UnionSetArray.jl")
include("LazyOperations/Rectification.jl")  # must come after UnionSet

# =======
# Aliases
# =======
include("Interfaces/aliases.jl")
include("Interfaces/AbstractArraySet.jl")

# =============================
# Conversions between set types
# =============================
include("Convert/convert.jl")

# ===========================
# Concrete operations on sets
# ===========================
include("ConcreteOperations/cartesian_product.jl")
include("ConcreteOperations/convex_hull.jl")
include("ConcreteOperations/difference.jl")
include("ConcreteOperations/distance.jl")
include("ConcreteOperations/exact_sum.jl")
include("ConcreteOperations/intersection.jl")
include("ConcreteOperations/isapprox.jl")
include("ConcreteOperations/isdisjoint.jl")
include("ConcreteOperations/isequal.jl")
include("ConcreteOperations/isequivalent.jl")
include("ConcreteOperations/isstrictsubset.jl")
include("ConcreteOperations/issubset.jl")
include("ConcreteOperations/linear_combination.jl")
include("ConcreteOperations/minkowski_difference.jl")
include("ConcreteOperations/minkowski_sum.jl")
include("Utils/samples.jl")

# =====================
# Approximations module
# =====================
include("Approximations/Approximations.jl")
# We export all symbols from Approximations.
# Note that the LazySets module is not supposed to depend on Approximations.
# It can, however, happen that we forget to add the `using` statements.
@reexport using .Approximations

# ==================================
# Plotting (requires Approximations)
# ==================================
include("Plotting/plot_recipes.jl")
include("Plotting/mesh.jl")

# ==========================
# Parallel-algorithms module
# ==========================
include("Parallel/Parallel.jl")

# ==========================
# Parallel-algorithms module
# ==========================
include("MatrixSets/MatrixSets.jl")

# ==============================
# Activate assertions by default
# ==============================
activate_assertions()

# ===================================================
# Load external packages on-demand (using 'Requires')
# ===================================================
include("init.jl")

end # module
