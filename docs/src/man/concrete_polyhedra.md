# Concrete Polyhedra

The focus of `LazySets.jl` is to wrap set representations and operations into
specialized types, delaying the evaluation of the result of an expression until
it is necessary.
However, sometimes it is desirable to perform an explicit computation.
For concrete operations with polyhedra we rely on the polyhedra manipulation
library [`Polyhedra.jl`](https://github.com/JuliaPolyhedra/Polyhedra.jl).

`Polyhedra.jl` provides a unified interface to well-known implementations of
polyhedral computations (which we also call *backends*) such as *CDD* or *LRS*.
See the complete list in [the documentation of
`Polyhedra.jl`](https://juliapolyhedra.github.io/Polyhedra.jl/latest/installation/#Getting-Libraries-1).
This is a great advantage because we can easily use a library that supports
floating point arithmetic, rational arithmetic, multiple precision, etc.
The libraries also include projection and elimination of variables through
Fourier-Motzkin.
If you are interested in specific numeric types different from the default
`Float64`, such as `Float32`, these may not be supported by the backend, in
which case Julia may automatically promote to, e.g., `Float64`.
As an example, *CDD* (which is used via the wrapper package
[`CDDLib.jl`](https://github.com/JuliaPolyhedra/CDDLib.jl)) can only be used
with numeric type `Float64` for floating-point arithmetic and with numeric type
`Rational{BigInt}` for exact arithmetic.

Below we give examples of operations that are performed using `Polyhedra.jl`.

```@contents
Pages = ["concrete_polyhedra.md"]
Depth = 3
```

## Creating polyhedra

To use the `Polyhedra.jl` interface, you need to load the package with `using Polyhedra`.
Let's create an H-representation object:

```@example concrete_polyhedra
using Plots, LazySets, Polyhedra, LinearAlgebra

A = [1. 1;1 -1;-1 0]
b = [1.,0,0]
H = Polyhedra.hrep(A, b)
```

It is used to instantiate a new polyhedron:

```@example concrete_polyhedra
p = polyhedron(H)
```

Now, `p` is of the generic type `Polyhedra.SimplePolyhedron{2,Float64, ...}`, where
`2` states for its ambient dimension, and `Float64` the numeric field. The remaining
fields specify the type of representation:

```@example concrete_polyhedra
typeof(p)
```

Observe that we can use a particular backend, such as the `CDD` library:

```@example concrete_polyhedra
using CDDLib

p = polyhedron(H, CDDLib.Library())
```

On the other hand, a `LazySets.HPolytope` object can be constructed from `p`:

```@example concrete_polyhedra
x = HPolytope(p)
x.constraints
```

Conversely, from a `HPolytope` we can build a polyhedron:

```@example concrete_polyhedra
y = polyhedron(x)
typeof(y)
```

Moreover, you can specify the backend with an extra argument.
For instance, we can use an exact representation through the
`Library(:exact)`:

```@example concrete_polyhedra
A, b = Rational{Int}[1 1;1 -1;-1 0], Rational{Int}[1,0,0]
p = HPolytope(A, b)

polyhedron(p; backend=CDDLib.Library(:exact))
```

## Methods

The utility methods available are convex hull, intersection and cartesian
product.
The dual representation as a list of vertices can be obtained with the
`vertices_list` function.

```@example concrete_polyhedra
p = HPolytope([LinearConstraint([1.0, 0.0], 1.0),
               LinearConstraint([0.0, 1.0], 1.0),
               LinearConstraint([-1.0, 0.0], 1.0),
               LinearConstraint([0.0, -1.0], 1.0)])

constraints_list(p)
```

```@example concrete_polyhedra
vertices_list(p)
```

For example, the concrete intersection of two polytopes is performed with the
`intersection` method.

```@example concrete_polyhedra
E = Ellipsoid(ones(2), Diagonal([2.0, 0.5]))
B = Ball1([2.5, 1.5], .8)

import LazySets.Approximations.overapproximate
polyoverapprox(x) = HPolytope(overapproximate(x, 1e-3).constraints)

Epoly = polyoverapprox(E)
Bpoly = polyoverapprox(B)
X = intersection(Epoly, Bpoly)

plot(E, 1e-3, aspectratio=1, alpha=0.4)
plot!(B, 1e-3, alpha=0.4)
plot!(X, 1e-3, alpha=0.4, color="black")
```

## Projections

Projection of high-dimensional polyhedra and elimination of variables can be
performed with the `eliminate` function, which supports three types of methods:
`:FourierMotzkin`, `:BlockElimination` and `:ProjectGenerators`.

For further details, see [the documentation of
`Polyhedra.jl`](https://juliapolyhedra.github.io/Polyhedra.jl/latest/projection/).
