# Interval Hulls

In this section we illustrate some interval hull operators
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

## Balls and singleton plots

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

To add plots to the same pair of axes we used the `plot!` functin. We have added some
points of the set which are farthest in some given directions.
Single points can be plotted using the `Singleton` type.
In the third line we plot the center of the ball with a cross marker.


```@example example_ih
plot!(Singleton(σ([1., 0], X)))
plot!(Singleton(σ([1., 1], X)))
plot!(Singleton(X.center), markershape=:x)
```

!!! note
    To see the list of available plot keyword arguments, use the `plotattr([attr])`
    function, where `attr` is the symbol `:Plot`, `:Series`, `:Axis` or `:Subplot`.

## Symmetric interval hull

The symmetric interval hull is defined in the `Approximations` module.

```@example example_ih
import LazySets.Approximations.symmetric_interval_hull

S = symmetric_interval_hull(X)
plot!(S, alpha=0.2)
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

## Convex hull overapproximation

Given two balls in the 2-norm, we build their convex hull and the overapproximation
using the symmetric interval hull:

```@example example_ih
X = Ball2(ones(2), 0.5)
Y = Ball2([-3,-.5], 0.8)
Z = CH(X, Y)

plot(X, 1e-3, aspectratio=1)
plot!(Y, 1e-3)
plot!(Z, 1e-3, alpha=0.2)
plot!(symmetric_interval_hull(Z), alpha=0.1)
```
