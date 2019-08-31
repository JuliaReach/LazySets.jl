__precompile__(true)

# main module for `LazySets.jl`
module LazySets

using Requires, SparseArrays, LinearAlgebra, Reexport, MathProgBase,
      GLPKMathProgInterface
using LinearAlgebra: checksquare
import LinearAlgebra: norm, ×
import Random
using Random: AbstractRNG, GLOBAL_RNG, SamplerType, shuffle
import InteractiveUtils: subtypes

export Arrays
export ×

# =======================
# Arrays auxiliary module
# =======================
include("Arrays/Arrays.jl")
using .Arrays

# ===================
# Auxiliary functions
# ===================
include("helper_functions.jl")
include("comparisons.jl")
include("macros.jl")
include("samples.jl")

# ==================
# Abstract set types
# ==================
include("Interfaces/LazySet.jl")
include("Interfaces/AbstractPolyhedron.jl")
include("HalfSpace.jl") # must be here to make LinearConstraint available
include("Interfaces/AbstractPolyhedron_functions.jl")
include("Interfaces/AbstractPolytope.jl")
include("Interfaces/AbstractCentrallySymmetric.jl")
include("Interfaces/AbstractCentrallySymmetricPolytope.jl")
include("Interfaces/AbstractZonotope.jl")
include("Interfaces/AbstractHyperrectangle.jl")
include("Interfaces/AbstractPolygon.jl")
include("Interfaces/AbstractSingleton.jl")
include("Interfaces/AbstractHPolygon.jl")

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

# ==================================
# Types representing non-convex sets
# ==================================
include("PolynomialZonotope.jl")

# =================================
# Types representing set operations
# =================================
include("LazyOperations/CartesianProduct.jl")
include("LazyOperations/Complement.jl")
include("LazyOperations/ConvexHull.jl")
include("LazyOperations/ExponentialMap.jl")
include("LazyOperations/Intersection.jl")
include("LazyOperations/LinearMap.jl")
include("LazyOperations/AffineMap.jl")  # must come after LinearMap
include("LazyOperations/MinkowskiSum.jl")
include("LazyOperations/ResetMap.jl")
include("LazyOperations/SymmetricIntervalHull.jl")
include("LazyOperations/Translation.jl")
include("LazyOperations/UnionSet.jl")
include("LazyOperations/Rectification.jl")  # must come after UnionSet

# =======
# Aliases
# =======
include("Interfaces/aliases.jl")

# =============================
# Conversions between set types
# =============================
include("convert.jl")

# =====================
# Approximations module
# =====================
include("Approximations/Approximations.jl")
# We export all symbols from Approximations.
# Note that the LazySets module is not supposed to depend on Approximations.
# It can, however, happen that we forget to add the `using` statements.
@reexport using .Approximations

# ===========================
# Concrete operations on sets
# ===========================
include("ConcreteOperations/convex_hull.jl")
include("ConcreteOperations/difference.jl")
include("ConcreteOperations/intersection.jl")
include("ConcreteOperations/isdisjoint.jl")
include("ConcreteOperations/issubset.jl")
include("ConcreteOperations/minkowski_difference.jl")
include("ConcreteOperations/minkowski_sum.jl")

# ========
# Plotting
# ========
include("Plotting/plot_recipes.jl")
include("Plotting/mesh.jl")

# ==========================
# Parallel algorithms module
# ==========================
include("Parallel/Parallel.jl")

# ===================================================
# Load external packages on-demand (using 'Requires')
# ===================================================
include("init.jl")

end # module
