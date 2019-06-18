__precompile__(true)

# main module for `LazySets.jl`
module LazySets

using Requires, SparseArrays, LinearAlgebra
using LinearAlgebra: checksquare
import LinearAlgebra: norm, ×
import Random
using Random: AbstractRNG, GLOBAL_RNG, SamplerType, shuffle
import InteractiveUtils: subtypes

export Approximations
export ×

# ===================
# Auxiliary functions
# ===================
include("helper_functions.jl")
include("comparisons.jl")
include("macros.jl")

# ==================
# Abstract set types
# ==================
include("LazySet.jl")
include("AbstractPolyhedron.jl")
include("HalfSpace.jl") # must be here to make LinearConstraint available
include("AbstractPolyhedron_functions.jl")
include("AbstractPolytope.jl")
include("AbstractCentrallySymmetric.jl")
include("AbstractCentrallySymmetricPolytope.jl")
include("AbstractZonotope.jl")
include("AbstractHyperrectangle.jl")
include("AbstractPolygon.jl")
include("AbstractSingleton.jl")
include("AbstractHPolygon.jl")

# =============================
# Types representing basic sets
# =============================
include("Ball1.jl")
include("Ball2.jl")
include("BallInf.jl")
include("Ballp.jl")
include("Ellipsoid.jl")
include("EmptySet.jl")
include("HPolygon.jl")
include("HPolygonOpt.jl")
include("HPolytope.jl")
include("HPolyhedron.jl")
include("Hyperplane.jl")
include("Hyperrectangle.jl")
include("Interval.jl")
include("Line.jl")
include("LineSegment.jl")
include("Singleton.jl")
include("Universe.jl")
include("VPolygon.jl")
include("VPolytope.jl")
include("ZeroSet.jl")
include("Zonotope.jl")
include("Complement.jl")

# ==================================
# Types representing non-convex sets
# ==================================
include("PolynomialZonotope.jl")

# ===================================================
# Algorithms to compute the convex hull of polygons
# ===================================================
include("concrete_convex_hull.jl")

# =================================
# Types representing set operations
# =================================
include("CartesianProduct.jl")
include("ConvexHull.jl")
include("ExponentialMap.jl")
include("Intersection.jl")
include("LinearMap.jl")
include("MinkowskiSum.jl")
include("ResetMap.jl")
include("SymmetricIntervalHull.jl")
include("Translation.jl")
include("UnionSet.jl")
include("Rectification.jl")

# =============================
# Conversions between set types
# =============================
include("intersection_helper.jl")
include("convert.jl")

# =====================
# Approximations module
# =====================
include("Approximations/Approximations.jl")

# ===========================
# Concrete operations on sets
# ===========================
include("concrete_intersection.jl")
include("is_intersection_empty.jl")
include("is_subset.jl")
include("difference.jl")

# =======
# Aliases
# =======
include("aliases.jl")

# ==========================
# Parallel algorithms module
# ==========================
include("Parallel/Parallel.jl")

# ============
# Plot recipes
# ============
include("plot_recipes.jl")

# ===================================================
# Load external packages on-demand (using 'Requires')
# ===================================================
include("init.jl")

end # module
