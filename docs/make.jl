using Documenter, LazySets

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
        "A Reachability Algorithm" => "man/reach_zonotopes.md"],
        "Library" => Any[
        "Set Interfaces" => "lib/interfaces.md",
        "Common Set Representations" => "lib/representations.md",
        "Common Set Operations" => "lib/operations.md",
        "Approximations" => "lib/approximations.md",
        "Utility Functions" => "lib/utils.md"],
        "About" => "about.md"
    ]
)

"deploy" in ARGS && deploydocs(
    repo = "github.com/JuliaReach/LazySets.jl.git",
    target = "build",
    julia  = "0.6",
    deps = nothing,
    make = nothing
)