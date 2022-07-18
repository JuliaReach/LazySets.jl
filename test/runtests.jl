using LazySets, LazySets.Approximations, Test, LinearAlgebra, SparseArrays,
      StaticArrays

# fix random number generator seed
using Random
Random.seed!(1234)

# ========================
# Optional dependencies
# ========================
import Distributions, ExponentialUtilities, Expokit, IntervalArithmetic,
       IntervalMatrices, Optim, Pkg, TaylorModels
const IA = IntervalArithmetic
using IntervalArithmetic: IntervalBox
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
global test_suite_doctests = VERSION >= v"1.7"
global test_suite_polyhedra = true
global test_suite_plotting = true

if (length(ARGS) == 0) || (ARGS[1] == "--default")
    # default test suite including doctests
elseif ARGS[1] == "--basic"
    # basic test suite
    test_suite_doctests = false
    test_suite_polyhedra = false
    test_suite_plotting = false
elseif ARGS[1] == "--polyhedra"
    # Polyhedra.jl test suite
    test_suite_doctests = false
    test_suite_polyhedra = true
    test_suite_plotting = false
elseif ARGS[1] == "--plot"
    # plotting test suite
    test_suite_doctests = false
    test_suite_polyhedra = false
    test_suite_plotting = true
elseif ARGS[1] == "--all"
    # complete test suite
    test_suite_doctests = true
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

    # =========================
    # Testing utility functions
    # =========================
    @time @testset "LazySets.Util" begin include("Utils/util.jl") end
    @time @testset "LazySets.Comparisons" begin include("Utils/comparisons.jl") end
    @time @testset "LazySets.Interior" begin include("Utils/interior.jl") end

    # ======================================================
    # Testing auxiliary machinery for vectors and matrices
    # ======================================================
    @time @testset "LazySets.Arrays" begin include("Arrays/array_operations.jl") end

    # =======================================
    # Testing interfaces to external packages
    # =======================================
    @time @testset "LazySets.CDDLib" begin include("Interfaces/CDDLib.jl") end

    # =======================================
    # Testing types that inherit from ConvexSet
    # =======================================
    @time @testset "LazySets.ConvexSet" begin include("Interfaces/ConvexSet.jl") end
    @time @testset "LazySets.Singleton" begin include("Sets/Singleton.jl") end
    @time @testset "LazySets.Ball1" begin include("Sets/Ball1.jl") end
    @time @testset "LazySets.Ball2" begin include("Sets/Ball2.jl") end
    @time @testset "LazySets.Ballp" begin include("Sets/Ballp.jl") end
    @time @testset "LazySets.BallInf" begin include("Sets/BallInf.jl") end
    @time @testset "LazySets.Hyperrectangle" begin include("Sets/Hyperrectangle.jl") end
    @time @testset "LazySets.Polygon" begin include("Sets/Polygon.jl") end
    @time @testset "LazySets.Polytope" begin include("Sets/Polytope.jl") end
    @time @testset "LazySets.Polyhedron" begin include("Sets/Polyhedron.jl") end
    @time @testset "LazySets.Zonotope" begin include("Sets/Zonotope.jl") end
    @time @testset "LazySets.ZeroSet" begin include("Sets/ZeroSet.jl") end
    @time @testset "LazySets.EmptySet" begin include("Sets/EmptySet.jl") end
    @time @testset "LazySets.Ellipsoid" begin include("Sets/Ellipsoid.jl") end
    @time @testset "LazySets.Hyperplane" begin include("Sets/Hyperplane.jl") end
    @time @testset "LazySets.HalfSpace" begin include("Sets/HalfSpace.jl") end
    @time @testset "LazySets.Interval" begin include("Sets/Interval.jl") end
    @time @testset "LazySets.LineSegment" begin include("Sets/LineSegment.jl") end
    @time @testset "LazySets.Line2D" begin include("Sets/Line2D.jl") end
    @time @testset "LazySets.Line" begin include("Sets/Line.jl") end
    @time @testset "LazySets.Universe" begin include("Sets/Universe.jl") end
    @time @testset "LazySets.HParallelotope" begin include("Sets/HParallelotope.jl") end
    @time @testset "LazySets.RotatedHyperrectangle" begin include("Sets/RotatedHyperrectangle.jl") end
    @time @testset "LazySets.Star" begin include("Sets/Star.jl") end

    # =========================================
    # Testing types representing set operations
    # =========================================
    @time @testset "LazySets.Bloating" begin include("LazyOperations/Bloating.jl") end
    @time @testset "LazySets.Intersection" begin include("LazyOperations/Intersection.jl") end
    @time @testset "LazySets.ConvexHull" begin include("LazyOperations/ConvexHull.jl") end
    @time @testset "LazySets.ExponentialMap" begin include("LazyOperations/ExponentialMap.jl") end
    @time @testset "LazySets.LinearMap" begin include("LazyOperations/LinearMap.jl") end
    @time @testset "LazySets.InverseLinearMap" begin include("LazyOperations/InverseLinearMap.jl") end
    @time @testset "LazySets.MinkowskiSum" begin include("LazyOperations/MinkowskiSum.jl") end
    @time @testset "LazySets.CartesianProduct" begin include("LazyOperations/CartesianProduct.jl") end
    @time @testset "LazySets.ResetMap" begin include("LazyOperations/ResetMap.jl") end
    @time @testset "LazySets.SymmetricIntervalHull" begin include("LazyOperations/SymmetricIntervalHull.jl") end
    @time @testset "LazySets.Translation" begin include("LazyOperations/Translation.jl") end
    @time @testset "LazySets.AffineMap" begin include("LazyOperations/AffineMap.jl") end

    # ======================
    # Testing set interfaces
    # ======================
    @time @testset "LazySets.CompactSet" begin include("Interfaces/CompactSet.jl") end

    # =========================================================
    # Testing other set types that do not inherit from ConvexSet
    # =========================================================
    @time @testset "LazySets.Complement" begin include("LazyOperations/Complement.jl") end
    @time @testset "LazySets.DensePolynomialZonotope" begin include("Sets/DensePolynomialZonotope.jl") end
    @time @testset "LazySets.Rectification" begin include("LazyOperations/Rectification.jl") end
    @time @testset "LazySets.SimpleSparsePolynomialZonotope" begin include("Sets/SimpleSparsePolynomialZonotope.jl") end
    @time @testset "LazySets.SparsePolynomialZonotope" begin include("Sets/SparsePolynomialZonotope.jl") end
    @time @testset "LazySets.UnionSet" begin include("LazyOperations/UnionSet.jl") end

    # ===================
    # Concrete operations
    # ===================
    @time @testset "LazySets.area" begin include("ConcreteOperations/area.jl") end
    @time @testset "LazySets.concrete_convex_hull" begin include("ConcreteOperations/convex_hull.jl") end
    @time @testset "LazySets.concrete_cartesian_product" begin include("ConcreteOperations/cartesian_product.jl") end
    @time @testset "LazySets.issubset" begin include("ConcreteOperations/issubset.jl") end
    @time @testset "LazySets.isdisjoint" begin include("ConcreteOperations/isdisjoint.jl") end
    @time @testset "LazySets.distance" begin include("ConcreteOperations/distance.jl") end
    @time @testset "LazySets.isstrictsubset" begin include("ConcreteOperations/isstrictsubset.jl") end
    @time @testset "LazySets.samples" begin include("ConcreteOperations/samples.jl") end

    # =================================================================
    # Algorithms for approximation of convex sets using support vectors
    # =================================================================
    @time @testset "LazySets.Approximations.overapproximation" begin include("Approximations/overapproximate.jl") end
    @time @testset "LazySets.Approximations.underapproximation" begin include("Approximations/underapproximate.jl") end
    @time @testset "LazySets.Approximations.template_directions" begin include("Approximations/template_directions.jl") end
    @time @testset "LazySets.Approximations.box_approximation" begin include("Approximations/box_approximation.jl") end
    @time @testset "LazySets.Approximations.ballinf_approximation" begin include("Approximations/ballinf_approximation.jl") end
    @time @testset "LazySets.Approximations.symmetric_interval_hull" begin include("Approximations/symmetric_interval_hull.jl") end
    @time @testset "LazySets.Approximations.radiusdiameter" begin include("Approximations/radiusdiameter.jl") end
    @time @testset "LazySets.Approximations.decompose" begin include("Approximations/decompose.jl") end
    @time @testset "LazySets.Approximations.hausdorff_distance" begin include("Approximations/hausdorff_distance.jl") end

    # ========================
    # Testing method ambiguity
    # ========================
    @time @testset "LazySets.method_ambiguities" begin
        for package in [LazySets, Approximations, Arrays, LazySets.Parallel]
            ambiguities = detect_ambiguities(package)
            @test isempty(ambiguities)
        end
    end

    # ====================================
    # Testing common API of all interfaces
    # (must be the last test because it
    #  loads Polyhedra.jl)
    # ====================================
    include("Utils/check_method_implementation.jl")
    @time @testset "LazySets.interfaces" begin include("Interfaces/interfaces.jl") end
end

if test_suite_plotting
    # define `plot` function as `RecipesBase.apply_recipe`
    import RecipesBase
    struct DummyBackend <: RecipesBase.AbstractBackend end
    struct DummyPlot <: RecipesBase.AbstractPlot{DummyBackend} end
    Base.length(::DummyPlot) = 0
    dict = Dict{Symbol, Any}(:plot_object => DummyPlot())
    plot(args...; kwargs...) = RecipesBase.apply_recipe(dict, args...; kwargs...)

    @time @testset "LazySets.plotting" begin include("Utils/plot.jl") end
end

if test_suite_doctests
    using Documenter
    include("../docs/init.jl")
    @time @testset "LazySets.doctests" begin doctest(LazySets) end
end
