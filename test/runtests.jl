using LazySets, LazySets.Approximations, Test, LinearAlgebra, SparseArrays

import IntervalArithmetic
const IA = IntervalArithmetic
using IntervalArithmetic: IntervalBox

# ========================
# Optional dependencies
# ========================
import Distributions, Expokit, IntervalMatrices, Optim, TaylorModels
using IntervalMatrices: ±, IntervalMatrix
using TaylorModels: set_variables, TaylorModelN

# ==============================
# Non-exported helper functions
# ==============================
using LazySets: ispermutation
using LazySets.Arrays: isinvertible, inner,
                       is_cyclic_permutation, SingleEntryVector

# conversion between numeric types
include("to_N.jl")

global test_suite_basic = true
global test_suite_doctests = true
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
    using LazySets: dim, HalfSpace, Interval, Line, translate
end

if test_suite_basic
    # =========================
    # Testing utility functions
    # =========================
    @time @testset "LazySets.Util" begin include("unit_util.jl") end
    @time @testset "LazySets.Comparisons" begin include("unit_comparisons.jl") end

    # =======================================
    # Testing types that inherit from LazySet
    # =======================================
    @time @testset "LazySets.Singleton" begin include("unit_Singleton.jl") end
    @time @testset "LazySets.Ball1" begin include("unit_Ball1.jl") end
    @time @testset "LazySets.Ball2" begin include("unit_Ball2.jl") end
    @time @testset "LazySets.Ballp" begin include("unit_Ballp.jl") end
    @time @testset "LazySets.BallInf" begin include("unit_BallInf.jl") end
    @time @testset "LazySets.Hyperrectangle" begin include("unit_Hyperrectangle.jl") end
    @time @testset "LazySets.Polygon" begin include("unit_Polygon.jl") end
    @time @testset "LazySets.Polytope" begin include("unit_Polytope.jl") end
    @time @testset "LazySets.Polyhedron" begin include("unit_Polyhedron.jl") end
    @time @testset "LazySets.Zonotope" begin include("unit_Zonotope.jl") end
    @time @testset "LazySets.ZeroSet" begin include("unit_ZeroSet.jl") end
    @time @testset "LazySets.EmptySet" begin include("unit_EmptySet.jl") end
    @time @testset "LazySets.Ellipsoid" begin include("unit_Ellipsoid.jl") end
    @time @testset "LazySets.Hyperplane" begin include("unit_Hyperplane.jl") end
    @time @testset "LazySets.HalfSpace" begin include("unit_HalfSpace.jl") end
    @time @testset "LazySets.Interval" begin include("unit_Interval.jl") end
    @time @testset "LazySets.LineSegment" begin include("unit_LineSegment.jl") end
    @time @testset "LazySets.Line" begin include("unit_Line.jl") end
    @time @testset "LazySets.Universe" begin include("unit_Universe.jl") end

    # =========================================
    # Testing types representing set operations
    # =========================================
    @time @testset "LazySets.Intersection" begin include("unit_Intersection.jl") end
    @time @testset "LazySets.ConvexHull" begin include("unit_ConvexHull.jl") end
    @time @testset "LazySets.ExponentialMap" begin include("unit_ExponentialMap.jl") end
    @time @testset "LazySets.LinearMap" begin include("unit_LinearMap.jl") end
    @time @testset "LazySets.MinkowskiSum" begin include("unit_MinkowskiSum.jl") end
    @time @testset "LazySets.CartesianProduct" begin include("unit_CartesianProduct.jl") end
    @time @testset "LazySets.ResetMap" begin include("unit_ResetMap.jl") end
    @time @testset "LazySets.SymmetricIntervalHull" begin include("unit_SymmetricIntervalHull.jl") end
    @time @testset "LazySets.concrete_convex_hull" begin include("unit_convex_hull.jl") end
    @time @testset "LazySets.AffineMap" begin include("unit_AffineMap.jl") end

    # ======================
    # Testing set interfaces
    # ======================
    @time @testset "LazySets.CompactSet" begin include("unit_CompactSet.jl") end

    # =========================================================
    # Testing other set types that do not inherit from LazySet
    # =========================================================
    @time @testset "LazySets.Complement" begin include("unit_Complement.jl") end
    @time @testset "LazySets.PolynomialZonotope" begin include("unit_PolynomialZonotope.jl") end
    @time @testset "LazySets.Rectification" begin include("unit_Rectification.jl") end
    @time @testset "LazySets.UnionSet" begin include("unit_UnionSet.jl") end

    # =================================================================
    # Algorithms for approximation of convex sets using support vectors
    # =================================================================
    @time @testset "LazySets.Approximations.overapproximation" begin include("unit_overapproximate.jl") end
    @time @testset "LazySets.Approximations.template_directions" begin include("unit_template_directions.jl") end
    @time @testset "LazySets.Approximations.box_approximation" begin include("unit_box_approximation.jl") end
    @time @testset "LazySets.Approximations.ballinf_approximation" begin include("unit_ballinf_approximation.jl") end
    @time @testset "LazySets.Approximations.radiusdiameter" begin include("unit_radiusdiameter.jl") end
    @time @testset "LazySets.Approximations.decompose" begin include("unit_decompose.jl") end

    # ========================
    # Testing method ambiguity
    # ========================
    include("check_method_ambiguity_binary.jl")
    @time @testset "LazySets.binary_operations" begin include("unit_binary_operations.jl") end

    # ====================================
    # Testing common API of all interfaces
    # (must be the last test because it
    #  loads Polyhedra.jl)
    # ====================================
    include("check_method_implementation.jl")
    @time @testset "LazySets.interfaces" begin include("unit_interfaces.jl") end
end

if test_suite_plotting
    import Plots
    using Plots: plot

    # fix namespace conflicts with Plots
    using LazySets: center

    @time @testset "LazySets.plotting" begin include("unit_plot.jl") end
end

if test_suite_doctests
    using Documenter
    include("../docs/init.jl")
    @time @testset "LazySets.doctests" begin doctest(LazySets) end
end
