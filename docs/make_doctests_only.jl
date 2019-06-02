using Documenter, LazySets
import Polyhedra, Optim, Expokit

makedocs(
    sitename = "LazySets.jl",
    modules = Module[LazySets, Approximations],
    source = "src/lib",
    doctest = true,
    strict = true
)
