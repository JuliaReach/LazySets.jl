using Documenter, LazySets

makedocs(
    modules = [LazySets],
    format = :html,
    sitename = "LazySets.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
        "Getting Started" => "man/getting_started.md",
        "Support Vectors" => "man/support_vectors.md",
        "Polyhedral Approximations" => "man/polyhedral_approximations.md"],
        "Library" => Any[
        "Common Set Representations" => "lib/representations.md", 
        "Common Set Operations" => "lib/operations.md",
        "Approximations" => "lib/approximations.md"],
        "About" => "about.md"
    ]
)

deploydocs(
    repo = "github.com/mforets/LazySets.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps = Deps.pip("mkdocs", "python-markdown-math"),
    make = nothing
)