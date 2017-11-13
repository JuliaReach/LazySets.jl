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
       Approximations

"""
    LazySet

Abstract type for a lazy set.

Every concrete `LazySet` must define a `σ(d, X)`, representing the support
vector of `X` in a given direction `d`, and `dim`, the ambient dimension
of the set `X`.
"""
abstract type LazySet end

# ============================
# Auxiliary types or functions
# ============================
include("LinearConstraints.jl")
include("helper_functions.jl")

# ===============================
# Types that inherit from LazySet
# ===============================
include("VoidSet.jl")
include("Singleton.jl")
include("Ball2.jl")
include("BallInf.jl")
include("Hyperrectangle.jl")
#include("Polyhedron.jl")  # optional (long startup time!)
include("HPolygon.jl")
include("HPolygonOpt.jl")
include("VPolygon.jl")

# =================================
# Types representing set operations
# =================================
include("ConvexHull.jl")
include("CartesianProduct.jl")
include("ExponentialMap.jl")
include("LinearMap.jl")
include("MinkowskiSum.jl")

# =================================================================
# Algorithms for approximation of convex sets using support vectors
# =================================================================
include("support_function.jl")
include("Approximations/Approximations.jl")
include("plot_recipes.jl")

end # module
