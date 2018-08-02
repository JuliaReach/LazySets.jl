using Documenter, LazySets

makedocs(
    doctest = true,
    modules = Module[LazySets, Approximations],
    source = "src/lib",
    strict = true
)
