using Test, SafeTestsets

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

@safetestset "API" begin
    include("API.jl")
end

using LazySets

# fix seed of random number generator (for reproducibility)
using Random
seed = max(abs(rand(Int)), 0)
println("using random seed $seed")
Random.seed!(seed)

# optional dependencies
@ts begin
    # only load the packages (no symbols)
    import CDDLib, Distributions, ExponentialUtilities, Expokit, GeometryBasics,
           IntervalBoxes, IntervalConstraintProgramming, IntervalMatrices,
           Ipopt, MiniQhull, Optim, PkgVersion, Polyhedra, RangeEnclosures, SCS,
           StaticArrays, TaylorModels
    if VERSION < v"1.12"
        import SetProg  # TODO add back unconditionally once it works in v1.12
    end

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

    @safetestset "CDDLib" begin
        include("Interfaces/CDDLib.jl")
    end
    @safetestset "SymEngine" begin
        include("Interfaces/SymEngine.jl")
    end

    # =================
    # LazySet interface
    # =================

    @safetestset "LazySet" begin
        include("Interfaces/LazySet.jl")
    end

    # ===============
    # Basic set types
    # ===============

    @safetestset "Ball1" begin
        include("Sets/Ball1.jl")
    end
    @safetestset "Ball2" begin
        include("Sets/Ball2.jl")
    end
    @safetestset "Ballp" begin
        include("Sets/Ballp.jl")
    end
    @safetestset "BallInf" begin
        include("Sets/BallInf.jl")
    end
    @safetestset "DensePolynomialZonotope" begin
        include("Sets/DensePolynomialZonotope.jl")
    end
    @safetestset "Ellipsoid" begin
        include("Sets/Ellipsoid.jl")
    end
    @safetestset "EmptySet" begin
        include("Sets/EmptySet.jl")
    end
    @safetestset "HalfSpace" begin
        include("Sets/HalfSpace.jl")
    end
    @safetestset "Hyperplane" begin
        include("Sets/Hyperplane.jl")
    end
    @safetestset "Hyperrectangle" begin
        include("Sets/Hyperrectangle.jl")
    end
    @safetestset "HParallelotope" begin
        include("Sets/HParallelotope.jl")
    end
    @safetestset "Interval" begin
        include("Sets/Interval.jl")
    end
    @safetestset "Line" begin
        include("Sets/Line.jl")
    end
    @safetestset "Line2D" begin
        include("Sets/Line2D.jl")
    end
    @safetestset "LineSegment" begin
        include("Sets/LineSegment.jl")
    end
    @safetestset "PolygonConvex" begin
        include("Sets/Polygon.jl")
    end
    @safetestset "PolygonNonconvex" begin
        include("Sets/PolygonNC.jl")
    end
    @safetestset "Polyhedron" begin
        include("Sets/Polyhedron.jl")
    end
    @safetestset "Polytope" begin
        include("Sets/Polytope.jl")
    end
    @safetestset "Universe" begin
        include("Sets/Universe.jl")
    end
    @safetestset "SimpleSparsePolynomialZonotope" begin
        include("Sets/SimpleSparsePolynomialZonotope.jl")
    end
    @safetestset "Singleton" begin
        include("Sets/Singleton.jl")
    end
    @safetestset "SparsePolynomialZonotope" begin
        include("Sets/SparsePolynomialZonotope.jl")
    end
    @safetestset "Star" begin
        include("Sets/Star.jl")
    end
    @safetestset "Tetrahedron" begin
        include("Sets/Tetrahedron.jl")
    end
    @safetestset "ZeroSet" begin
        include("Sets/ZeroSet.jl")
    end
    @safetestset "Zonotope" begin
        include("Sets/Zonotope.jl")
    end
    @safetestset "ZonotopeMD" begin
        include("Sets/ZonotopeMD.jl")
    end

    # =====================
    # Matrix set operations
    # =====================

    @safetestset "MatrixZonotope" begin
        include("MatrixSets/MatrixZonotope.jl")
    end

    # ===================
    # Lazy set operations
    # ===================

    @safetestset "AffineMap" begin
        include("LazyOperations/AffineMap.jl")
    end
    @safetestset "Bloating" begin
        include("LazyOperations/Bloating.jl")
    end
    @safetestset "CartesianProduct" begin
        include("LazyOperations/CartesianProduct.jl")
    end
    @safetestset "Complement" begin
        include("LazyOperations/Complement.jl")
    end
    @safetestset "ConvexHull" begin
        include("LazyOperations/ConvexHull.jl")
    end
    @safetestset "ExponentialMap" begin
        include("LazyOperations/ExponentialMap.jl")
    end
    @safetestset "ExactSum" begin
        include("LazyOperations/ExactSum.jl")
    end
    @safetestset "Intersection" begin
        include("LazyOperations/Intersection.jl")
    end
    @safetestset "InverseLinearMap" begin
        include("LazyOperations/InverseLinearMap.jl")
    end
    @safetestset "LinearMap" begin
        include("LazyOperations/LinearMap.jl")
    end
    @safetestset "MinkowskiSum" begin
        include("LazyOperations/MinkowskiSum.jl")
    end
    @safetestset "Rectification" begin
        include("LazyOperations/Rectification.jl")
    end
    @safetestset "ResetMap" begin
        include("LazyOperations/ResetMap.jl")
    end
    @safetestset "SymmetricIntervalHull" begin
        include("LazyOperations/SymmetricIntervalHull.jl")
    end
    @safetestset "Translation" begin
        include("LazyOperations/Translation.jl")
    end
    @safetestset "UnionSet" begin
        include("LazyOperations/UnionSet.jl")
    end

    # ==============
    # Set interfaces
    # ==============

    @safetestset "CompactSet" begin
        include("Interfaces/CompactSet.jl")
    end
    @safetestset "AbstractZonotope" begin
        include("Interfaces/AbstractZonotope.jl")
    end

    # =======================
    # Concrete set operations
    # =======================

    @safetestset "area" begin
        include("ConcreteOperations/area.jl")
    end
    @safetestset "cartesian_product" begin
        include("ConcreteOperations/cartesian_product.jl")
    end
    @safetestset "convex_hull" begin
        include("ConcreteOperations/convex_hull.jl")
    end
    @safetestset "difference" begin
        include("ConcreteOperations/difference.jl")
    end
    @safetestset "distance" begin
        include("ConcreteOperations/distance.jl")
    end
    @safetestset "exact_sum" begin
        include("ConcreteOperations/exact_sum.jl")
    end
    @safetestset "interior" begin
        include("ConcreteOperations/interior.jl")
    end
    @safetestset "intersection" begin
        include("ConcreteOperations/intersection.jl")
    end
    @safetestset "isdisjoint" begin
        include("ConcreteOperations/isdisjoint.jl")
    end
    @safetestset "isequivalent" begin
        include("ConcreteOperations/isequivalent.jl")
    end
    @safetestset "isstrictsubset" begin
        include("ConcreteOperations/isstrictsubset.jl")
    end
    @safetestset "issubset" begin
        include("ConcreteOperations/issubset.jl")
    end
    @safetestset "linear_combination" begin
        include("ConcreteOperations/linear_combination.jl")
    end
    @safetestset "minkowski_difference" begin
        include("ConcreteOperations/minkowski_difference.jl")
    end
    @safetestset "minkowski_sum" begin
        include("ConcreteOperations/minkowski_sum.jl")
    end
    @safetestset "samples" begin
        include("ConcreteOperations/samples.jl")
    end

    # =====================
    # Approximations module
    # =====================

    @safetestset "ballinf_approximation" begin
        include("Approximations/ballinf_approximation.jl")
    end
    @safetestset "box_approximation" begin
        include("Approximations/box_approximation.jl")
    end
    @safetestset "decompose" begin
        include("Approximations/decompose.jl")
    end
    @safetestset "hausdorff_distance" begin
        include("Approximations/hausdorff_distance.jl")
    end
    @safetestset "overapproximate" begin
        include("Approximations/overapproximate.jl")
    end
    @safetestset "overapproximate_norm" begin
        include("Approximations/overapproximate_norm.jl")
    end
    @safetestset "overapproximate_expmap" begin
        include("Approximations/overapproximate_expmap.jl")
    end
    @safetestset "overapproximate_matrixzonotope" begin
        include("Approximations/overapproximate_matrixzonotope.jl")
    end
    @safetestset "radiusdiameter" begin
        include("Approximations/radiusdiameter.jl")
    end
    @safetestset "symmetric_interval_hull" begin
        include("Approximations/symmetric_interval_hull.jl")
    end
    @safetestset "template_directions" begin
        include("Approximations/template_directions.jl")
    end
    @safetestset "underapproximate" begin
        include("Approximations/underapproximate.jl")
    end

    # ====================
    # Solver functionality
    # ====================

    @safetestset "lp_solvers" begin
        include("Utils/lp_solvers.jl")
    end

    # ============================
    # Common API of all interfaces
    # ============================

    @safetestset "interfaces" begin
        include("Interfaces/interfaces.jl")
    end

    # ========
    # Plotting
    # ========

    @safetestset "plotting" begin
        include("Utils/plot.jl")
    end

    # =================
    # Quality assurance
    # =================

    include("quality_assurance.jl")
end
