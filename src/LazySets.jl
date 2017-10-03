__precompile__(true)
"""
Main module for `LazySets.jl` -- a Julia package for calculus with convex sets.


Every LazySet type must define:
 σ(d, sf) -- Vector{Float64}, support vector in a given direction
 dim      -- Int64, dimension
"""
module LazySets

abstract type LazySet end

# ============================
# Auxiliary types or functions
# ============================
include("LinearConstraints.jl")
include("HelperFuncs.jl")

# ===============================
# Types that inherit from LazySet
# ===============================
include("VoidSet.jl")
include("Singleton.jl")
include("Ball2.jl")
include("BallInf.jl")
include("Hyperrectangle.jl")
#include("LazySets/Polyhedron.jl")  # optional
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

export LazySet, dim, σ, support_vector, ρ, support_function

# =================================================================
# Algorithms for approximation of convex sets using support vectors
# =================================================================
include("Approximations.jl")

end
