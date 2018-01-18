__precompile__(true)

"""
Main module for `LazySets.jl` -- a Julia package for calculus with convex sets.
"""
module LazySets

using RecipesBase, IterTools

export LazySet,
       ρ, support_function,
       σ, support_vector,
       dim,
       norm,
       radius,
       diameter,
       Approximations

"""
    LazySet{N}

Abstract type for a lazy set.

### Notes

Every concrete `LazySet` must define the following functions:
- `σ(d::AbstractVector{N}, S::LazySet)::AbstractVector{N}` -- the support vector
    of `S` in a given direction `d`
- `dim(S::LazySet)::Int` -- the ambient dimension of `S`

`LazySet` types should be parameterized with a type `N`, typically `N<:Real`,
for using different numeric types.
"""
abstract type LazySet{N} end

# ============================
# Auxiliary types or functions
# ============================
include("LinearConstraints.jl")
include("helper_functions.jl")

# ===============================
# Types that inherit from LazySet
# ===============================
# abstract types
include("AbstractPolytope.jl")
include("AbstractPointSymmetric.jl")
include("AbstractPointSymmetricPolytope.jl")
include("AbstractHyperrectangle.jl")
include("AbstractPolygon.jl")
include("AbstractHPolygon.jl")
include("AbstractSingleton.jl")
# concrete types
include("EmptySet.jl")
include("ZeroSet.jl")
include("Singleton.jl")
include("Ball2.jl")
include("BallInf.jl")
include("Ball1.jl")
include("Ballp.jl")
include("Hyperrectangle.jl")
include("HPolytope.jl")
include("HPolygon.jl")
include("HPolygonOpt.jl")
include("VPolygon.jl")
include("Zonotope.jl")
include("Ellipsoid.jl")

# =================================
# Types representing set operations
# =================================
include("ConvexHull.jl")
include("CartesianProduct.jl")
include("ExponentialMap.jl")
include("LinearMap.jl")
include("MinkowskiSum.jl")
include("Intersection.jl")
include("SymmetricIntervalHull.jl")

# =================================================================
# Algorithms for approximation of convex sets using support vectors
# =================================================================
include("support_function.jl")
include("Approximations/Approximations.jl")
include("plot_recipes.jl")

end # module
