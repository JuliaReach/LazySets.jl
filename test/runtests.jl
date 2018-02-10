#!/usr/bin/env julia
using LazySets, Base.Test

# conversion between numeric types
include("to_N.jl")

# =======================================
# Testing types that inherit from LazySet
# =======================================
@time @testset "LazySets.Singleton" begin include("unit_Singleton.jl") end
@time @testset "LazySets.Ball2" begin include("unit_Ball2.jl") end
@time @testset "LazySets.BallInf" begin include("unit_BallInf.jl") end
@time @testset "LazySets.Hyperrectangle" begin include("unit_Hyperrectangle.jl") end
@time @testset "LazySets.Polygon" begin include("unit_Polygon.jl") end
@time @testset "LazySets.Polytope" begin include("unit_Polytope.jl") end
@time @testset "LazySets.Zonotope" begin include("unit_Zonotope.jl") end
@time @testset "LazySets.ZeroSet" begin include("unit_ZeroSet.jl") end
@time @testset "LazySets.EmptySet" begin include("unit_EmptySet.jl") end
@time @testset "LazySets.Ellipsoid" begin include("unit_Ellipsoid.jl") end
@time @testset "LazySets.Hyperplane" begin include("unit_Hyperplane.jl") end
@time @testset "LazySets.HalfSpace" begin include("unit_HalfSpace.jl") end

# =========================================
# Testing types representing set operations
# =========================================
@time @testset "LazySets.Intersection" begin include("unit_Intersection.jl") end
@time @testset "LazySets.ConvexHull" begin include("unit_ConvexHull.jl") end
@time @testset "LazySets.ExponentialMap" begin include("unit_ExponentialMap.jl") end
@time @testset "LazySets.LinearMap" begin include("unit_LinearMap.jl") end
@time @testset "LazySets.MinkowskiSum" begin include("unit_MinkowskiSum.jl") end
@time @testset "LazySets.CartesianProduct" begin include("unit_CartesianProduct.jl") end
@time @testset "LazySets.SymmetricIntervalHull" begin include("unit_SymmetricIntervalHull.jl") end

# =================================================================
# Algorithms for approximation of convex sets using support vectors
# =================================================================
@time @testset "LazySets.Approximations.overapproximation" begin include("unit_overapproximation2D.jl") end
@time @testset "LazySets.Approximations.box_approximation" begin include("unit_box_approximation.jl") end
@time @testset "LazySets.Approximations.ballinf_approximation" begin include("unit_ballinf_approximation.jl") end
@time @testset "LazySets.Approximations.radiusdiameter" begin include("unit_radiusdiameter.jl") end
@time @testset "LazySets.Approximations.decompose" begin include("unit_decompose.jl") end

# ==================
# Further algorithms
# ==================
@time @testset "LazySets.binary_operations" begin include("unit_binary_operations.jl") end

# ====================================
# Testing common API of all interfaces
# (must be the last test because it
#  loads Polyhedra.jl)
# ====================================
@time @testset "LazySets.interfaces" begin include("unit_interfaces.jl") end
