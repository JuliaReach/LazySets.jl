using Documenter, LazySets

makedocs(
    #doctest = true,  # uncomment locally for skipping doctests (saves time!)
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
        "Iterative refinement" => "man/iterative_refinement.md"],
        "Library" => Any[
        "Common Set Representations" => "lib/representations.md", 
        "Common Set Operations" => "lib/operations.md",
        "Approximations" => "lib/approximations.md",
        "Utility Functions" => "lib/utils.md"],
        "About" => "about.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaReach/LazySets.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps = Deps.pip("mkdocs", "python-markdown-math"),
    make = nothing
)