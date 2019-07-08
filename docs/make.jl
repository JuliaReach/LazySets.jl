ENV["GKSwstype"] = "100"  # set 'GR environment' to 'no output' (for Travis CI)
using Documenter, LazySets
import Polyhedra, Optim, Expokit, TaylorModels

makedocs(
    sitename = "LazySets.jl",
    modules = Module[LazySets, Approximations],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/juliareach.css"]),
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
        "Getting Started" => "man/getting_started.md",
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
        "Set Interfaces" => "lib/interfaces.md",
        "Common Set Representations" => "lib/representations.md",
        "Common Set Operations" => "lib/operations.md",
        "Comparisons" => "lib/comparisons.md",
        "Conversions between set representations" => "lib/conversion.md",
        "Binary Functions on Sets" => "lib/binary_functions.md",
        "Approximations" => "lib/approximations.md",
        "Utility Functions" => "lib/utils.md",
        "Parallel" => "lib/parallel.md"
        ],
        "About" => "about.md"
    ],
    doctest = false,
    strict = true
)

deploydocs(
    repo = "github.com/JuliaReach/LazySets.jl.git"
)
