using Documenter, LazySets, Polyhedra, Optim

makedocs(
    doctest = true,
    modules = Module[LazySets, Approximations],
    source = "src/lib",
    strict = true
)
