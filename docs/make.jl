using Documenter, LazySets
import Polyhedra, Optim, ExponentialUtilities, TaylorModels, Distributions,
       MiniQhull, Symbolics

include("init.jl")

makedocs(
    sitename = "LazySets.jl",
    modules = [LazySets, Approximations, Arrays, LazySets.Parallel],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
                             assets = ["assets/aligned.css"]),
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "Getting Started" => "man/getting_started.md",
            "Optional Features" => "man/optional_dependencies.md",
            "A Tour of LazySets" => "man/tour.md",
            "Introduction to Convex Sets" => "man/convex_sets.md",
            "Polyhedral Approximations" => "man/polyhedral_approximations.md",
            "Decomposing an Affine Map" => "man/decompose_example.md",
            "Fast 2D LPs" => "man/fast_2d_LPs.md",
            "Iterative Refinement" => "man/iterative_refinement.md",
            "Interval Hulls" => "man/interval_hulls.md",
            "Convex Hulls" => "man/convex_hulls.md",
            "Operations on Sets" => "man/set_operations.md",
            "A Reachability Algorithm" => "man/reach_zonotopes.md",
            "A Hybrid Reachability Algorithm" => "man/reach_zonotopes_hybrid.md",
            "Concrete Polyhedra" => "man/concrete_polyhedra.md",
            "Parallel Approximations" => "man/parallel_approximations.md",
            "Lazy Intersections" => "man/lazy_intersections.md"
        ],
        "Library" => Any[
            "Interfaces" => "lib/interfaces.md",
            "Sets" => [
                "Ball1" => "lib/sets/Ball1.md",
                "Ball2" => "lib/sets/Ball2.md",
                "BallInf" => "lib/sets/BallInf.md",
                "Ballp" => "lib/sets/Ballp.md",
                "Ellipsoid" => "lib/sets/Ellipsoid.md",
                "EmptySet" => "lib/sets/EmptySet.md",
                "HalfSpace" => "lib/sets/HalfSpace.md",
                "HParallelotope" => "lib/sets/HParallelotope.md",
                "HPolygon" => "lib/sets/HPolygon.md",
                "HPolygonOpt" => "lib/sets/HPolygonOpt.md",
                "HPolyhedron" => "lib/sets/HPolyhedron.md",
                "HPolytope" => "lib/sets/HPolytope.md",
                "Hyperplane" => "lib/sets/Hyperplane.md",
                "Hyperrectangle" => "lib/sets/Hyperrectangle.md",
                "Interval" => "lib/sets/Interval.md",
                "Line2D" => "lib/sets/Line2D.md",
                "Line" => "lib/sets/Line.md",
                "LineSegment" => "lib/sets/LineSegment.md",
                "DensePolynomialZonotope" => "lib/sets/DensePolynomialZonotope.md",
                "RotatedHyperrectangle" => "lib/sets/RotatedHyperrectangle.md",
                "SimpleSparsePolynomialZonotope" => "lib/sets/SimpleSparsePolynomialZonotope.md",
                "Singleton" => "lib/sets/Singleton.md",
                "Star" => "lib/sets/Star.md",
                "Universe" => "lib/sets/Universe.md",
                "VPolygon" => "lib/sets/VPolygon.md",
                "VPolytope" => "lib/sets/VPolytope.md",
                "ZeroSet" => "lib/sets/ZeroSet.md",
                "Zonotope" => "lib/sets/Zonotope.md",
            ],
            "Lazy Operations" => [
                "AffineMap" => "lib/lazy_operations/AffineMap.md",
                "Bloating" => "lib/lazy_operations/Bloating.md",
                "CartesianProduct" => "lib/lazy_operations/CartesianProduct.md",
                "Complement" => "lib/lazy_operations/Complement.md",
                "ConvexHull" => "lib/lazy_operations/ConvexHull.md",
                "ExponentialMap" => "lib/lazy_operations/ExponentialMap.md",
                "Intersection" => "lib/lazy_operations/Intersection.md",
                "LinearMap" => "lib/lazy_operations/LinearMap.md",
                "InverseLinearMap" => "lib/lazy_operations/InverseLinearMap.md",
                "MinkowskiSum" => "lib/lazy_operations/MinkowskiSum.md",
                "QuadraticMap" => "lib/lazy_operations/QuadraticMap.md",
                "Rectification" => "lib/lazy_operations/Rectification.md",
                "ResetMap" => "lib/lazy_operations/ResetMap.md",
                "SymmetricIntervalHull" => "lib/lazy_operations/SymmetricIntervalHull.md",
                "Translation" => "lib/lazy_operations/Translation.md",
                "UnionSet" => "lib/lazy_operations/UnionSet.md",
            ],
            "Concrete Operations" => "lib/binary_functions.md",
            "Conversions between set representations" => "lib/conversion.md",
            "Comparisons" => "lib/comparisons.md",
            "Approximations" => "lib/approximations.md",
            "Utility Functions" => "lib/utils.md",
            "Parallel" => "lib/parallel.md",
        ],
        "About" => "about.md"
    ],
    doctest = false,
    strict = true
)

deploydocs(
    repo = "github.com/JuliaReach/LazySets.jl.git",
    push_preview=true
)
