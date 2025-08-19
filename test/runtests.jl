using Test

# auxiliary code to skip expensive tests
begin
    __test_short = haskey(ENV, "JULIA_PKGEVAL")

    __test_Float64_only = false

    macro ts(arg)
        if !__test_short
            return quote
                $(esc(arg))
            end
        end
    end

    macro tv(v1, v2)
        if __test_short
            return v1
        else
            return @eval vcat($v1, $v2)
        end
    end

    macro tN(v)
        if __test_Float64_only
            return quote
                [$(esc(v))[1]]
            end
        else
            return v
        end
    end
end

@testset "LazySets.API" begin
    include("API.jl")
end

using LazySets

# fix seed of random number generator (for reproducibility)
using Random
Random.seed!(1234)

# optional dependencies
@ts begin
    # only load the packages (no symbols)
    import CDDLib, Distributions, ExponentialUtilities, Expokit, IntervalConstraintProgramming,
           IntervalMatrices, Ipopt, MiniQhull, Optim, PkgVersion, Polyhedra, RangeEnclosures, SCS,
           SetProg, StaticArrays, TaylorModels

    # load all symbols from the following packages
    using SymEngine, Symbolics
    # resolve conflict
    using LazySets: factors
end

# tests
@ts begin
    # ==============================
    # Interface to external packages
    # ==============================

    @testset "LazySets.CDDLib" begin
        include("Interfaces/CDDLib.jl")
    end
    @testset "LazySets.SymEngine" begin
        include("Interfaces/SymEngine.jl")
    end

    # =================
    # LazySet interface
    # =================

    @testset "LazySets.LazySet" begin
        include("Interfaces/LazySet.jl")
    end

    # ===============
    # Basic set types
    # ===============

    @testset "LazySets.Ball1" begin
        include("Sets/Ball1.jl")
    end
    @testset "LazySets.Ball2" begin
        include("Sets/Ball2.jl")
    end
    @testset "LazySets.Ballp" begin
        include("Sets/Ballp.jl")
    end
    @testset "LazySets.BallInf" begin
        include("Sets/BallInf.jl")
    end
    @testset "LazySets.DensePolynomialZonotope" begin
        include("Sets/DensePolynomialZonotope.jl")
    end
    @testset "LazySets.Ellipsoid" begin
        include("Sets/Ellipsoid.jl")
    end
    @testset "LazySets.EmptySet" begin
        include("Sets/EmptySet.jl")
    end
    @testset "LazySets.HalfSpace" begin
        include("Sets/HalfSpace.jl")
    end
    @testset "LazySets.Hyperplane" begin
        include("Sets/Hyperplane.jl")
    end
    @testset "LazySets.Hyperrectangle" begin
        include("Sets/Hyperrectangle.jl")
    end
    @testset "LazySets.HParallelotope" begin
        include("Sets/HParallelotope.jl")
    end
    @testset "LazySets.Interval" begin
        include("Sets/Interval.jl")
    end
    @testset "LazySets.Line" begin
        include("Sets/Line.jl")
    end
    @testset "LazySets.Line2D" begin
        include("Sets/Line2D.jl")
    end
    @testset "LazySets.LineSegment" begin
        include("Sets/LineSegment.jl")
    end
    @testset "LazySets.PolygonConvex" begin
        include("Sets/Polygon.jl")
    end
    @testset "LazySets.PolygonNonconvex" begin
        include("Sets/PolygonNC.jl")
    end
    @testset "LazySets.Polyhedron" begin
        include("Sets/Polyhedron.jl")
    end
    @testset "LazySets.Polytope" begin
        include("Sets/Polytope.jl")
    end
    @testset "LazySets.Universe" begin
        include("Sets/Universe.jl")
    end
    @testset "LazySets.SimpleSparsePolynomialZonotope" begin
        include("Sets/SimpleSparsePolynomialZonotope.jl")
    end
    @testset "LazySets.Singleton" begin
        include("Sets/Singleton.jl")
    end
    @testset "LazySets.SparsePolynomialZonotope" begin
        include("Sets/SparsePolynomialZonotope.jl")
    end
    @testset "LazySets.Star" begin
        include("Sets/Star.jl")
    end
    @testset "LazySets.Tetrahedron" begin
        include("Sets/Tetrahedron.jl")
    end
    @testset "LazySets.ZeroSet" begin
        include("Sets/ZeroSet.jl")
    end
    @testset "LazySets.Zonotope" begin
        include("Sets/Zonotope.jl")
    end
    @testset "LazySets.ZonotopeMD" begin
        include("Sets/ZonotopeMD.jl")
    end

    # =====================
    # Matrix set operations
    # =====================

    @testset "LazySets.MatrixZonotope" begin
        include("MatrixSets/MatrixZonotope.jl")
    end

    # ===================
    # Lazy set operations
    # ===================

    @testset "LazySets.AffineMap" begin
        include("LazyOperations/AffineMap.jl")
    end
    @testset "LazySets.Bloating" begin
        include("LazyOperations/Bloating.jl")
    end
    @testset "LazySets.CartesianProduct" begin
        include("LazyOperations/CartesianProduct.jl")
    end
    @testset "LazySets.Complement" begin
        include("LazyOperations/Complement.jl")
    end
    @testset "LazySets.ConvexHull" begin
        include("LazyOperations/ConvexHull.jl")
    end
    @testset "LazySets.ExponentialMap" begin
        include("LazyOperations/ExponentialMap.jl")
    end
    @testset "LazySets.Intersection" begin
        include("LazyOperations/Intersection.jl")
    end
    @testset "LazySets.InverseLinearMap" begin
        include("LazyOperations/InverseLinearMap.jl")
    end
    @testset "LazySets.LinearMap" begin
        include("LazyOperations/LinearMap.jl")
    end
    @testset "LazySets.MinkowskiSum" begin
        include("LazyOperations/MinkowskiSum.jl")
    end
    @testset "LazySets.Rectification" begin
        include("LazyOperations/Rectification.jl")
    end
    @testset "LazySets.ResetMap" begin
        include("LazyOperations/ResetMap.jl")
    end
    @testset "LazySets.SymmetricIntervalHull" begin
        include("LazyOperations/SymmetricIntervalHull.jl")
    end
    @testset "LazySets.Translation" begin
        include("LazyOperations/Translation.jl")
    end
    @testset "LazySets.UnionSet" begin
        include("LazyOperations/UnionSet.jl")
    end

    # ==============
    # Set interfaces
    # ==============

    @testset "LazySets.CompactSet" begin
        include("Interfaces/CompactSet.jl")
    end
    @testset "LazySets.AbstractZonotope" begin
        include("Interfaces/AbstractZonotope.jl")
    end

    # =======================
    # Concrete set operations
    # =======================

    @testset "LazySets.area" begin
        include("ConcreteOperations/area.jl")
    end
    @testset "LazySets.cartesian_product" begin
        include("ConcreteOperations/cartesian_product.jl")
    end
    @testset "LazySets.convex_hull" begin
        include("ConcreteOperations/convex_hull.jl")
    end
    @testset "LazySets.difference" begin
        include("ConcreteOperations/difference.jl")
    end
    @testset "LazySets.distance" begin
        include("ConcreteOperations/distance.jl")
    end
    @testset "LazySets.exact_sum" begin
        include("ConcreteOperations/exact_sum.jl")
    end
    @testset "LazySets.interior" begin
        include("ConcreteOperations/interior.jl")
    end
    @testset "LazySets.intersection" begin
        include("ConcreteOperations/intersection.jl")
    end
    @testset "LazySets.isdisjoint" begin
        include("ConcreteOperations/isdisjoint.jl")
    end
    @testset "LazySets.isequivalent" begin
        include("ConcreteOperations/isequivalent.jl")
    end
    @testset "LazySets.isstrictsubset" begin
        include("ConcreteOperations/isstrictsubset.jl")
    end
    @testset "LazySets.issubset" begin
        include("ConcreteOperations/issubset.jl")
    end
    @testset "LazySets.minkowski_difference" begin
        include("ConcreteOperations/minkowski_difference.jl")
    end
    @testset "LazySets.minkowski_sum" begin
        include("ConcreteOperations/minkowski_sum.jl")
    end
    @testset "LazySets.samples" begin
        include("ConcreteOperations/samples.jl")
    end

    # =====================
    # Approximations module
    # =====================

    @testset "LazySets.Approximations.ballinf_approximation" begin
        include("Approximations/ballinf_approximation.jl")
    end
    @testset "LazySets.Approximations.box_approximation" begin
        include("Approximations/box_approximation.jl")
    end
    @testset "LazySets.Approximations.decompose" begin
        include("Approximations/decompose.jl")
    end
    @testset "LazySets.Approximations.hausdorff_distance" begin
        include("Approximations/hausdorff_distance.jl")
    end
    @testset "LazySets.Approximations.overapproximate" begin
        include("Approximations/overapproximate.jl")
    end
    @testset "LazySets.Approximations.overapproximate_norm" begin
        include("Approximations/overapproximate_norm.jl")
    end
    @testset "LazySets.Approximations.overapproximate_expmap" begin
        include("Approximations/overapproximate_expmap.jl")
    end
    @testset "LazySets.Approximations.radiusdiameter" begin
        include("Approximations/radiusdiameter.jl")
    end
    @testset "LazySets.Approximations.symmetric_interval_hull" begin
        include("Approximations/symmetric_interval_hull.jl")
    end
    @testset "LazySets.Approximations.template_directions" begin
        include("Approximations/template_directions.jl")
    end
    @testset "LazySets.Approximations.underapproximate" begin
        include("Approximations/underapproximate.jl")
    end

    # ====================
    # Solver functionality
    # ====================

    @testset "LazySets.lp_solvers" begin
        include("Utils/lp_solvers.jl")
    end

    # ============================
    # Common API of all interfaces
    # ============================

    include("Utils/check_method_implementation.jl")
    @testset "LazySets.interfaces" begin
        include("Interfaces/interfaces.jl")
    end

    # ========
    # Plotting
    # ========

    @testset "LazySets.plotting" begin
        include("Utils/plot.jl")
    end

    # ====
    # Aqua
    # ====

    include("Aqua.jl")
end
