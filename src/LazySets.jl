module LazySets

using Reexport

include("API/API.jl")
@reexport using .API
import .API: eltype, extrema, isdisjoint, isempty, \, ∈, ≈, ==, ⊆,
             rand, norm, permute, distance, rectify,
             affine_map, an_element, area, center, complement, concretize,
             constraints_list, constraints, convex_hull, diameter, dim, exponential_map,
             high, is_interior_point, is_polyhedral, isbounded, isboundedtype,
             isconvextype, isempty, isoperation, isoperationtype, isuniversal,
             linear_map, low, norm, project, radius, reflect, sample, scale,
             scale!, support_function, ρ, support_vector, σ, surface, translate,
             translate!, vertices_list, vertices, volume,
             cartesian_product, difference, distance, exact_sum, intersection,
             is_intersection_empty, isequivalent, ⊂, linear_combination,
             minkowski_difference, pontryagin_difference, minkowski_sum

import Base: copy, rationalize, \
import LinearAlgebra: ×, normalize, normalize!
import RecipesBase: apply_recipe

export Arrays
export ×, normalize, subtypes

using LinearAlgebra, RecipesBase, Requires, SparseArrays
import GLPK, JuMP, Random, ReachabilityBase
import IntervalArithmetic as IA
using LinearAlgebra: checksquare
using Random: AbstractRNG, GLOBAL_RNG, SamplerType, shuffle, randperm

@static if VERSION < v"1.9"
    stack(vecs) = hcat(vecs...)
end

# ================
# ReachabilityBase
# ================

import ReachabilityBase.Assertions
using ReachabilityBase.Assertions: @assert
include("Utils/assertions.jl")

using ReachabilityBase.Require
using ReachabilityBase.Comparison
using ReachabilityBase.Iteration
using ReachabilityBase.Commutative
using ReachabilityBase.Distribution
using ReachabilityBase.Subtypes
using ReachabilityBase.Arrays
using ReachabilityBase.Basetype

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
include("Interfaces/AbstractPolyhedron.jl")
include("Sets/HalfSpace.jl")  # must come before AbstractPolyhedron_functions
include("Interfaces/AbstractPolyhedron_functions.jl")
include("Interfaces/AbstractPolytope.jl")
include("Interfaces/AbstractCentrallySymmetric.jl")
include("Interfaces/AbstractCentrallySymmetricPolytope.jl")
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
include("Sets/Universe.jl")

include("Sets/EmptySet/EmptySetModule.jl")
@reexport using ..EmptySetModule: EmptySet, ∅, _isdisjoint_emptyset

include("Sets/Ball1.jl")
include("Sets/Ball2.jl")
include("Sets/BallInf.jl")
include("Sets/Ballp.jl")
include("Sets/DensePolynomialZonotope.jl")
include("Sets/Ellipsoid.jl")
include("Sets/HParallelotope.jl")
include("Sets/HPolygon.jl")
include("Sets/HPolygonOpt.jl")
include("Sets/HPolytope.jl")
include("Sets/HPolyhedron.jl")
include("Sets/Hyperplane.jl")
include("Sets/Hyperrectangle.jl")
include("Sets/Line2D.jl")
include("Sets/Line.jl")
include("Sets/LineSegment.jl")
include("Sets/Polygon.jl")
include("Sets/RotatedHyperrectangle.jl")
include("Sets/Singleton.jl")
include("Sets/SimpleSparsePolynomialZonotope.jl")
include("Sets/SparsePolynomialZonotope.jl")
include("Sets/VPolygon.jl")
include("Sets/VPolytope.jl")
include("Sets/Tetrahedron.jl")
include("Sets/ZeroSet.jl")
include("Sets/Zonotope.jl")

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
include("Sets/Star.jl")

# =============================
# Conversions between set types
# =============================
include("convert.jl")

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
