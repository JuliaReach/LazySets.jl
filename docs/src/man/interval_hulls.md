# Interval Hulls

In this section we illustrate the interval hull operators
as well as several plotting functionalities.

```@contents
Pages = ["interval_hulls.md"]
Depth = 3
```

```@meta
DocTestSetup = quote
    using LazySets, Plots, LazySets.Approximations
end
```

## Balls and Singletons

Consider a ball in the 2-norm. By default, the coefficients of this set are 64-bit floating point
numbers. Other numeric types (such as lower precision floating point, or rational)
can be defined with the proper argument types in the `Ball2` constructor.

```@example example_ih
using LazySets, Plots

X = Ball2(ones(2), 0.5)
```

To plot a lazy set, we use the `plot` function. Recall that, by design, lazy sets
plots overapproximate with box directions only. To have a sharp definition of the
borders, use the accuracy as a second argument.

```@example example_ih
plot(X, 1e-3, aspectratio=1)
```

To add plots to the same pair of axes we use `plot!`. Let's add some points of
the set which are farthest in some given directions. Single points can be plotted
using the `Singleton` type.
In the third line we plot the center of the ball picking a custom cross marker.


```@example example_ih
plot!(Singleton(σ([1., 0], X)))
plot!(Singleton(σ([1., 1], X)))
plot!(Singleton(X.center), markershape=:x)
```

!!! note
    To see the list of available plot keyword arguments, use the `plotattr([attr])`
    function, where `attr` is the symbol `:Plot`, `:Series`, `:Axis` or `:Subplot`.

For the remaining of this section we define another ball in the 2-norm and their
convex hull.

```@example example_ih
Y = Ball2([-3,-.5], 0.8)
Z = CH(X, Y)

plot(X, 1e-3, aspectratio=1)
plot!(Y, 1e-3)
plot!(Z, 1e-3, alpha=0.2)
```

## Ballinf approximation

A simple overapproximation with a `BallInf` is obtained with the `ballinf_approximation`
function. It overapproximates a convex set by a tight ball in the infinity norm by
evaluating the support vector in the canonical directions.

```@example example_ih
import LazySets.Approximations.ballinf_approximation

plot(X, 1e-3, aspectratio=1)
plot!(Y, 1e-3)
plot!(Z, 1e-3, alpha=0.2)

Bapprox = ballinf_approximation(Z)

plot!(Bapprox, alpha=0.1)
plot!(Singleton(Bapprox.center), markershape=:x)
```

```@example example_ih
Bapprox.center, Bapprox.radius
```

## Interval hull approximation

If we want to have different lengths for each dimension, instead of the
`ballinf_approximation`, we can use the approximation with a hyperrectangle through
the `interval_hull` function.

```@example example_ih
import LazySets.Approximations.interval_hull

plot(X, 1e-3, aspectratio=1)
plot!(Y, 1e-3)
plot!(Z, 1e-3, alpha=0.2)

Happrox = interval_hull(Z)

plot!(Happrox, alpha=0.1)
plot!(Singleton(Happrox.center), markershape=:x)
```

```@example example_ih
Happrox.center, Happrox.radius
```

!!! note
    The `interval_hull` function is an alias for the `box_approximation` function.
    The nomenclature for approximation functions is `*_approximation_*`. To see a
    list of all approximation functions, either search in the docs or type
    `names(LazySets.Approximations)`.

## Symmetric interval hull

Contrary to the previous approximations, the symmetric interval hull is centered
around the origin. It is defined in the `Approximations` module as well.

```@example example_ih
import LazySets.Approximations.symmetric_interval_hull

plot(X, 1e-3, aspectratio=1)
plot!(Y, 1e-3)
plot!(Z, 1e-3, alpha=0.2)

S = symmetric_interval_hull(Z)
plot!(S, alpha=0.2)
plot!(Singleton(S.center), markershape=:x)
```

```@example example_ih
S.center, S.radius
```

We can get the list of vertices using the `vertices_list` function:

```@example example_ih
vertices_list(S)
```
 
For instance, compute the support vector in the south-east direction:
 
```@example example_ih
σ([1., -1.], S)
```
 
 It is also possible to pass a sparse vector as direction, and the result is a
 sparse vector:
 
 ```@example example_ih
σ(sparsevec([1., -1.]), S)
```

## Norm, radius and diameter

In this part we illustrate some functions to obtain metric properties of sets.
These functions apply generally to any `LazySet`. For some types, specialized
methods are triggered automatically through multiple-dispatch.

The *norm* of a convex set is the norm of the enclosing ball (of the given norm)
of minimal volume. For instance:

```@example example_ih
import LazySets.Approximations: norm, radius, diameter

norm(X), norm(Y), norm(Z)
```

The *radius* of a convex set. It is the radius of the enclosing ball
(of the given norm) of minimal volume with the same center. In the previous example,

```@example example_ih
radius(X), radius(Y), radius(Z)
```

Finally, it is sometimes convenient to ask directly the diameter of the set,
defined as twice the radius:

```@example example_ih
diameter(X), diameter(Y), diameter(Z)
```

In the previous examples, we have used the infinity norm (default). It is possible
to pass the required `p`-norm as a second argument.
