__precompile__(true)

"""
Main module for `LazySets.jl` -- a Julia package for calculus with convex sets.
"""
module LazySets

using RecipesBase, IterTools, Requires

export Approximations

# ============================
# Auxiliary types or functions
# ============================
include("LinearConstraints.jl")
include("helper_functions.jl")

# ==================
# Abstract set types
# ==================
include("LazySet.jl")
include("AbstractPolytope.jl")
include("AbstractPointSymmetric.jl")
include("AbstractPointSymmetricPolytope.jl")
include("AbstractHyperrectangle.jl")
include("AbstractPolygon.jl")
include("AbstractHPolygon.jl")
include("AbstractSingleton.jl")

# =============================
# Types representing basic sets
# =============================
include("Ball1.jl")
include("Ball2.jl")
include("BallInf.jl")
include("Ballp.jl")
include("Ellipsoid.jl")
include("EmptySet.jl")
include("HalfSpace.jl")
include("HPolygon.jl")
include("HPolygonOpt.jl")
include("HPolytope.jl")
include("Hyperplane.jl")
include("Hyperrectangle.jl")
include("Singleton.jl")
include("VPolygon.jl")
include("VPolytope.jl")
include("ZeroSet.jl")
include("Zonotope.jl")
include("PolynomialZonotope.jl")

# =================================
# Types representing set operations
# =================================
include("CartesianProduct.jl")
include("ConvexHull.jl")
include("ExponentialMap.jl")
include("Intersection.jl")
include("LinearMap.jl")
include("MinkowskiSum.jl")
include("SymmetricIntervalHull.jl")

# =============================
# Conversions between set types
# =============================
include("convert.jl")

# =======================================
# Algorithms for check operations on sets
# =======================================
include("is_intersection_empty.jl")
include("is_subset.jl")

# =====================
# Approximations module
# =====================
include("Approximations/Approximations.jl")

# ============
# Plot recipes
# ============
include("plot_recipes.jl")

end # module
