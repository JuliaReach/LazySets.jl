ENV["GKSwstype"] = "100"  # set 'GR environment' to 'no output' (for Travis CI)
using Documenter, LazySets, Polyhedra

makedocs(
    doctest = true,  # use this flag to skip doctests (saves time!)
    modules = [LazySets, Approximations],
    format = :html,
    assets = ["assets/juliareach.css"],
    sitename = "LazySets.jl",
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
        ],
        "Library" => Any[
        "Set Interfaces" => "lib/interfaces.md",
        "Common Set Representations" => "lib/representations.md",
        "Common Set Operations" => "lib/operations.md",
        "Conversion between set representations" => "lib/conversion.md",
        "Binary Functions on Sets" => "lib/binary_functions.md",
        "Approximations" => "lib/approximations.md",
        "Utility Functions" => "lib/utils.md",
#         "Methods Collection" => "lib/methods_fix.md",
        ],
        "About" => "about.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaReach/LazySets.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps = nothing,
    make = nothing
)
