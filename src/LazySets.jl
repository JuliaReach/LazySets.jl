__precompile__(true)

"""
Main module for `LazySets.jl` -- a Julia package for calculus with convex sets.
"""
module LazySets

export LazySet, dim, σ, support_vector, ρ, support_function,
       Approximations

"""
    LazySet

Abstract type for a lazy set. Every concrete `LazySet` must define:

    ``σ(d, X)`` -- support vector of ``X`` in a given direction ``d``
      dim       -- the ambient dimension
"""
abstract type LazySet end

# ============================
# Auxiliary types or functions
# ============================
include("LinearConstraintaints.jl")
include("HelperFuncs.jl")

# ===============================
# Types that inherit from LazySet
# ===============================
include("VoidSet.jl")
include("Singleton.jl")
include("Ball2.jl")
include("BallInf.jl")
include("Hyperrectangle.jl")
#include("Polyhedron.jl")  # optional (long time)
include("Polygon.jl")

# =================================
# Types representing set operations
# =================================
include("ConvexHull.jl")
include("CartesianProduct.jl")
include("ExponentialMap.jl")
include("LinearMap.jl")
include("MinkowskiSum.jl")

"""
    ρ(d::Vector{Float64}, sf::LazySet)::Float64

Evaluate the support function of a set in a given direction.

# INPUT:

- `d`  -- a real vector, the direction investigated
- `sf` -- a convex set

# OUTPUT:

- `ρ(d, sf)` -- the support function
"""
function ρ(d::Vector{Float64}, sf::LazySet)::Float64
    return dot(d, σ(d, sf))
end

# alias
support_function = ρ
support_vector = σ

# =================================================================
# Algorithms for approximation of convex sets using support vectors
# =================================================================
include("Approximations.jl")

end # module
