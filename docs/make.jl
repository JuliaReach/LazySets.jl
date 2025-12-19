ENV["GKSwstype"] = "100"  # prevent plots from opening interactively

using Documenter, LazySets, DocumenterCitations
import Plots, Polyhedra, Optim, ExponentialUtilities, TaylorModels, Distributions,
       MiniQhull, Symbolics, SymEngine, IntervalMatrices

include("init.jl")

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:alpha)

makedocs(; sitename="LazySets.jl",
         modules=[LazySets, LazySets.API, Approximations, LazySets.Parallel],
         format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true",
                                assets=["assets/aligned.css", "assets/citations.css"],
                                size_threshold_warn=150 * 2^10),
         pagesonly=true,
         plugins=[bib],
         pages=[
                #
                "Home" => "index.md",
                "Manual" => [
                             #
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
                             "Unary Operations on Sets" => "man/unary_set_operations.md",
                             "A Reachability Algorithm" => "man/reach_zonotopes.md",
                             "A Hybrid Reachability Algorithm" => "man/reach_zonotopes_hybrid.md",
                             "Concrete Polyhedra" => "man/concrete_polyhedra.md",
                             "Parallel Approximations" => "man/parallel_approximations.md",
                             "Lazy Intersections" => "man/lazy_intersections.md"
                             #
                             ],
                "Library" => ["API" => "lib/API.md",
                              "Set Interfaces" => [
                                                   #
                                                   "lib/interfaces/overview.md",
                                                   "lib/interfaces/LazySet.md",
                                                   "lib/interfaces/ConvexSet.md",
                                                   "lib/interfaces/AbstractCentrallySymmetric.md",
                                                   "lib/interfaces/AbstractPolyhedron.md",
                                                   "lib/interfaces/AbstractPolytope.md",
                                                   "lib/interfaces/AbstractPolygon.md",
                                                   "lib/interfaces/AbstractHPolygon.md",
                                                   "lib/interfaces/AbstractCentrallySymmetricPolytope.md",
                                                   "lib/interfaces/AbstractZonotope.md",
                                                   "lib/interfaces/AbstractHyperrectangle.md",
                                                   "lib/interfaces/AbstractSingleton.md",
                                                   "lib/interfaces/AbstractAffineMap.md",
                                                   "lib/interfaces/AbstractPolynomialZonotope.md",
                                                   "lib/interfaces/AbstractSparsePolynomialZonotope.md",
                                                   "lib/interfaces/AbstractBallp.md"
                                                   #
                                                   ],
                              "Sets" => [
                                         #
                                         "Ball1" => "lib/sets/Ball1.md",
                                         "Ball2" => "lib/sets/Ball2.md",
                                         "BallInf" => "lib/sets/BallInf.md",
                                         "Ballp" => "lib/sets/Ballp.md",
                                         "DensePolynomialZonotope" => "lib/sets/DensePolynomialZonotope.md",
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
                                         "Polygon" => "lib/sets/Polygon.md",
                                         "SimpleSparsePolynomialZonotope" => "lib/sets/SimpleSparsePolynomialZonotope.md",
                                         "SparsePolynomialZonotope" => "lib/sets/SparsePolynomialZonotope.md",
                                         "Singleton" => "lib/sets/Singleton.md",
                                         "Star" => "lib/sets/Star.md",
                                         "Tetrahedron" => "lib/sets/Tetrahedron.md",
                                         "Universe" => "lib/sets/Universe.md",
                                         "VPolygon" => "lib/sets/VPolygon.md",
                                         "VPolytope" => "lib/sets/VPolytope.md",
                                         "ZeroSet" => "lib/sets/ZeroSet.md",
                                         "Zonotope" => "lib/sets/Zonotope.md",
                                         "ZonotopeMD" => "lib/sets/ZonotopeMD.md"
                                         #
                                         ],
                              "Lazy Operations" => [
                                                    #
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
                                                    "UnionSet" => "lib/lazy_operations/UnionSet.md"
                                                    #
                                                    ],
                              "Concrete Binary Operations" => [
                                                               #
                                                               "lib/concrete_binary_operations/overview.md",
                                                               "lib/concrete_binary_operations/cartesian_product.md",
                                                               "lib/concrete_binary_operations/convex_hull.md",
                                                               "lib/concrete_binary_operations/difference.md",
                                                               "lib/concrete_binary_operations/distance.md",
                                                               "lib/concrete_binary_operations/intersection.md",
                                                               "lib/concrete_binary_operations/minkowski_sum.md",
                                                               "lib/concrete_binary_operations/minkowski_difference.md",
                                                               "lib/concrete_binary_operations/isdisjoint.md",
                                                               "lib/concrete_binary_operations/issubset.md"
                                                               #
                                                               ],
                              "Conversions between set representations" => "lib/conversion.md",
                              "Approximations" => [
                                                   #
                                                   "lib/approximations/overview.md",
                                                   "lib/approximations/overapproximate.md",
                                                   "lib/approximations/box_approximation.md",
                                                   "lib/approximations/iterative_refinement.md",
                                                   "lib/approximations/template_directions.md",
                                                   "lib/approximations/underapproximate.md",
                                                   "lib/approximations/approximate.md",
                                                   "lib/approximations/decompose.md",
                                                   "lib/approximations/hausdorff_distance.md",
                                                   "lib/approximations/overapproximate_norm.md",
                                                   "lib/approximations/overapproximate_expmap.md"
                                                   #
                                                   ],
                              "MatrixSets" => [
                                               #
                                               "MatrixZonotope" => "lib/matrixsets/MatrixZonotope.md"
                                               #
                                               ],
                              "Utilities" => "lib/utils.md",
                              "Parallel" => "lib/parallel.md"
                              #
                              ],
                "Bibliography" => "bibliography.md",
                "About" => "about.md"
                #
                ])

deploydocs(; repo="github.com/JuliaReach/LazySets.jl.git",
           push_preview=true)
