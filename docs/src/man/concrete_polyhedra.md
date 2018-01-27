# Concrete Polyhedra

The focus of `LazySets.jl` is to wrap set representations and operations into
specialized types, delaying the evaluation of the result of an expression until it
is necessary. However, sometimes it is necessary to do an explicit computation.
For concrete operations with polyhedra we rely on the polyhedra manipulation library
[Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl).

Actually, `Polyhedra.jl` provides a unified interface to well-known
implementations of polyhedral computations, such as CDD, PPL or LRS (see the
complete list [in the documentation of `Polyhedra.jl`](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation.html#Getting-Libraries-1)).
This is a great advantage because we can easily use a library that support floating
point arithmetic, rational arithmetic, multiple precision, etc.
The libraries also include projection and elimination of variables through Fourier-Motzkin.

Below we give examples of operations that are actually done via `Polyhedra.jl`.

```@contents
Pages = ["concrete_polyhedra.md"]
Depth = 3
```

```@meta
DocTestSetup = quote
    using LazySets, Plots, LazySets.Approximations, Polyhedra
end
```

## Creating polyhedra

We create a 2D polytope with the `HPolytope` constructor:

```@example concrete_polyhedra
using LazySets, Plots, Polyhedra

A = [1. 1;1 -1;-1 0]
b = [1.,0,0]
hrep = SimpleHRepresentation(A, b)
p = polyhedron(hrep)
```

The `HPolytope` type can be constructed from a `HRep` polyhedron:

```@example concrete_polyhedra
x = HPolytope(p)
x.constraints
```

Conversely, from a `HPolytope` we can build a `HRep` polyhedron:

```@example concrete_polyhedra
y = polyhedron(x)
typeof(y)
```

## Methods

For example, the intersection of two polytopes is performed with the `intersect`
method.

```@example concrete_polyhedra
E = Ellipsoid(ones(2), diagm([2.0, 0.5]))
B = Ball1([2.5, 1.5], .8)

import LazySets.Approximations.overapproximate
polyoverapprox(x) = HPolytope(overapproximate(x, 1e-3).constraints_list)

Epoly = polyoverapprox(E)
Bpoly = polyoverapprox(B)
X = intersect(Epoly, Bpoly)

plot(E, 1e-3, aspectratio=1, alpha=0.4)
plot!(B, 1e-3, alpha=0.4)
plot!(X, 1e-3, alpha=0.4, color="black")
```
