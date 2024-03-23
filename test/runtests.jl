using LazySets, LazySets.Approximations, Test, LinearAlgebra, SparseArrays,
      StaticArrays, GLPK
using JuMP: optimizer_with_attributes

# fix random number generator seed
using Random
Random.seed!(1234)

# ========================
# Optional dependencies
# ========================
import Distributions, ExponentialUtilities, Expokit, IntervalMatrices, Ipopt,
       MiniQhull, Optim, RangeEnclosures, SCS, SetProg, TaylorModels
import IntervalArithmetic as IA
using IntervalArithmetic: IntervalBox, interval
@static if VERSION >= v"1.9"
    vIA = pkgversion(IA)
    vGLPK = pkgversion(GLPK)
else
    import PkgVersion
    vIA = PkgVersion.Version(IA)
    vGLPK = PkgVersion.Version(GLPK)
end
using IntervalMatrices: ±, IntervalMatrix
using TaylorModels: set_variables, TaylorModelN
# ICP currently leads to unsatisfiable package requirements
# using IntervalConstraintProgramming
using Symbolics

# ==============================
# Non-exported helper functions
# ==============================
using LazySets: _leq, _geq, isapproxzero, _isapprox, _ztol, ispermutation
using LazySets.Arrays: isinvertible, inner, allequal,
                       is_cyclic_permutation, SingleEntryVector

global test_suite_basic = true
global test_suite_polyhedra = true
global test_suite_plotting = true

if (length(ARGS) == 0) || (ARGS[1] == "--default")
    # default test suite
elseif ARGS[1] == "--basic"
    # basic test suite
    test_suite_polyhedra = false
    test_suite_plotting = false
elseif ARGS[1] == "--polyhedra"
    # Polyhedra.jl test suite
    test_suite_polyhedra = true
    test_suite_plotting = false
elseif ARGS[1] == "--plot"
    # plotting test suite
    test_suite_polyhedra = false
    test_suite_plotting = true
elseif ARGS[1] == "--all"
    # complete test suite
    test_suite_polyhedra = true
    test_suite_plotting = true
else
    error("unknown parameter 1: $(ARGS[1])")
end

if test_suite_polyhedra || test_suite_plotting
    import Polyhedra
    using CDDLib # for tests that require CDDLib specific backend=...

    # fix namespace conflicts with Polyhedra
    using LazySets: dim, HalfSpace, Interval, Line2D, translate
end

if test_suite_basic
    # =======================================
    # Testing interfaces to external packages
    # =======================================
    @testset "LazySets.CDDLib" begin
        include("Interfaces/CDDLib.jl")
    end

    # =======================
    # Testing basic set types
    # =======================
    @testset "LazySets.Singleton" begin
        include("Sets/Singleton.jl")
    end
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
    @testset "LazySets.Hyperrectangle" begin
        include("Sets/Hyperrectangle.jl")
    end
    @testset "LazySets.PolygonConvex" begin
        include("Sets/Polygon.jl")
    end
    @testset "LazySets.PolygonNonconvex" begin
        include("Sets/PolygonNC.jl")
    end
    @testset "LazySets.Polytope" begin
        include("Sets/Polytope.jl")
    end
    @testset "LazySets.Tetrahedron" begin
        include("Sets/Tetrahedron.jl")
    end
    @testset "LazySets.Polyhedron" begin
        include("Sets/Polyhedron.jl")
    end
    @testset "LazySets.Zonotope" begin
        include("Sets/Zonotope.jl")
    end
    @testset "LazySets.ZeroSet" begin
        include("Sets/ZeroSet.jl")
    end
    @testset "LazySets.EmptySet" begin
        include("Sets/EmptySet.jl")
    end
    @testset "LazySets.Ellipsoid" begin
        include("Sets/Ellipsoid.jl")
    end
    @testset "LazySets.Hyperplane" begin
        include("Sets/Hyperplane.jl")
    end
    @testset "LazySets.HalfSpace" begin
        include("Sets/HalfSpace.jl")
    end
    @testset "LazySets.Interval" begin
        include("Sets/Interval.jl")
    end
    @testset "LazySets.LineSegment" begin
        include("Sets/LineSegment.jl")
    end
    @testset "LazySets.Line2D" begin
        include("Sets/Line2D.jl")
    end
    @testset "LazySets.Line" begin
        include("Sets/Line.jl")
    end
    @testset "LazySets.Universe" begin
        include("Sets/Universe.jl")
    end
    @testset "LazySets.HParallelotope" begin
        include("Sets/HParallelotope.jl")
    end
    @testset "LazySets.RotatedHyperrectangle" begin
        include("Sets/RotatedHyperrectangle.jl")
    end
    @testset "LazySets.Star" begin
        include("Sets/Star.jl")
    end
    @testset "LazySets.DensePolynomialZonotope" begin
        include("Sets/DensePolynomialZonotope.jl")
    end
    @testset "LazySets.SimpleSparsePolynomialZonotope" begin
        include("Sets/SimpleSparsePolynomialZonotope.jl")
    end
    @testset "LazySets.SparsePolynomialZonotope" begin
        include("Sets/SparsePolynomialZonotope.jl")
    end

    # =========================================
    # Testing types representing set operations
    # =========================================
    @testset "LazySets.Bloating" begin
        include("LazyOperations/Bloating.jl")
    end
    @testset "LazySets.Intersection" begin
        include("LazyOperations/Intersection.jl")
    end
    @testset "LazySets.ConvexHull" begin
        include("LazyOperations/ConvexHull.jl")
    end
    @testset "LazySets.ExponentialMap" begin
        include("LazyOperations/ExponentialMap.jl")
    end
    @testset "LazySets.LinearMap" begin
        include("LazyOperations/LinearMap.jl")
    end
    @testset "LazySets.InverseLinearMap" begin
        include("LazyOperations/InverseLinearMap.jl")
    end
    @testset "LazySets.MinkowskiSum" begin
        include("LazyOperations/MinkowskiSum.jl")
    end
    @testset "LazySets.CartesianProduct" begin
        include("LazyOperations/CartesianProduct.jl")
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
    @testset "LazySets.AffineMap" begin
        include("LazyOperations/AffineMap.jl")
    end
    @testset "LazySets.Complement" begin
        include("LazyOperations/Complement.jl")
    end
    @testset "LazySets.Rectification" begin
        include("LazyOperations/Rectification.jl")
    end
    @testset "LazySets.UnionSet" begin
        include("LazyOperations/UnionSet.jl")
    end

    # ======================
    # Testing set interfaces
    # ======================
    @testset "LazySets.CompactSet" begin
        include("Interfaces/CompactSet.jl")
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
    @testset "LazySets.Interior" begin
        include("ConcreteOperations/interior.jl")
    end
    @testset "LazySets.intersection" begin
        include("ConcreteOperations/isstrictsubset.jl")
    end
    @testset "LazySets.isdisjoint" begin
        include("ConcreteOperations/isdisjoint.jl")
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

    # ====================================
    # Algorithms for approximation of sets
    # ====================================
    @testset "LazySets.Approximations.overapproximation" begin
        include("Approximations/overapproximate.jl")
    end
    @testset "LazySets.Approximations.underapproximation" begin
        include("Approximations/underapproximate.jl")
    end
    @testset "LazySets.Approximations.template_directions" begin
        include("Approximations/template_directions.jl")
    end
    @testset "LazySets.Approximations.box_approximation" begin
        include("Approximations/box_approximation.jl")
    end
    @testset "LazySets.Approximations.ballinf_approximation" begin
        include("Approximations/ballinf_approximation.jl")
    end
    @testset "LazySets.Approximations.symmetric_interval_hull" begin
        include("Approximations/symmetric_interval_hull.jl")
    end
    @testset "LazySets.Approximations.radiusdiameter" begin
        include("Approximations/radiusdiameter.jl")
    end
    @testset "LazySets.Approximations.decompose" begin
        include("Approximations/decompose.jl")
    end
    @testset "LazySets.Approximations.hausdorff_distance" begin
        include("Approximations/hausdorff_distance.jl")
    end

    # ================================
    # Testing shared utility functions
    # ================================
    @testset "LazySets.lp_solvers" begin
        include("Utils/lp_solvers.jl")
    end

    # =====================================================
    # Testing common API of all interfaces
    # (must be the last test because it loads Polyhedra.jl)
    # =====================================================
    include("Utils/check_method_implementation.jl")
    @testset "LazySets.interfaces" begin
        include("Interfaces/interfaces.jl")
    end
end

if test_suite_plotting
    # define `plot` function as `RecipesBase.apply_recipe`
    import RecipesBase
    struct DummyBackend <: RecipesBase.AbstractBackend end
    struct DummyPlot <: RecipesBase.AbstractPlot{DummyBackend} end
    Base.length(::DummyPlot) = 0
    dict = Dict{Symbol,Any}(:plot_object => DummyPlot())
    plot(args...; kwargs...) = RecipesBase.apply_recipe(dict, args...; kwargs...)

    @testset "LazySets.plotting" begin
        include("Utils/plot.jl")
    end
end

include("Aqua.jl")
