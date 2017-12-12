# Interval Hulls

In this section we illustrate some interval hull operators
as well as several plotting functionalities.

```@meta
DocTestSetup = quote
    using LazySets, Plots, LazySets.Approximations
end
```

## Balls and singleton plots

Consider a ball in the 2-norm:

```@example
X = Ball2(ones(2), 0.5)
```
Observe that, by default, the coefficients of this set are 64-bit floating point
numbers.

We can plot this set using the `plot` function. Recall that, by design, lazy sets
plots overapproximate with box directions only. To have a sharp definition of the
borders, use the accuracy as a second argument:

```julia
plot(X, 1e-3, aspectratio=1)
```

Let's add to this plot some points of the set which are farthest in some given directions.
Single points can be plotted using the `Singleton` type.

To add a new plot to the previous one, use `plot!`. In the third line we plot the
center of the ball with a cross marker.

!!! note
    To see the list of available plot keyword arguments, use the `plotattr([attr])`
    function, where `attr` is the symbol `:Plot`, `:Series`, `:Axis` or `:Subplot`.

```julia
plot!(Singleton(σ([1., 0], X)))
plot!(Singleton(σ([1., 1], X)))
plot!(Singleton(X.center), markershape=:x)
```

## Symmetric interval hull

The symmetric interval hull is defined in the `Approximations` module:

```julia
julia> S = Approximations.symmetric_interval_hull(X)
LazySets.Hyperrectangle{Float64}([0.0, 0.0], [1.5, 1.5])
```

```julia
plot!(S, alpha=0.2)
```

We can get the list of vertices using the `vertices_list` function:

```julia
julia> vertices_list(S)
4-element Array{Array{Float64,1},1}:
 [1.5, 1.5]
 [-1.5, 1.5]
 [1.5, -1.5]
 [-1.5, -1.5]
```
 
For instance, compute the support vector in the south-east direction:
 
```julia
σ([1., -1.], S)
2-element Array{Float64,1}:
  1.5
 -1.5
 ```
 
 It is also possible to pass a sparse vector as direction, and the result is a
 sparse vector:
 
 ```julia
julia> σ(sparsevec([1., -1.]), S)
2-element SparseVector{Float64,Int64} with 2 stored entries:
  [1]  =  1.5
  [2]  =  -1.5
```

## Convex hull overapproximation

Given two balls in the 2-norm, we build their convex hull and the overapproximation
using the symmetric interval hull:

```julia
X = Ball2(ones(2), 0.5)
Y = Ball2([-3,-.5], 0.8)
Z = CH(X, Y)

plot(X, 1e-3, aspectratio=1)
plot!(Y, 1e-3)
plot!(Z, 1e-3, alpha=0.2)

plot!(Approximations.symmetric_interval_hull(Z), alpha=0.1)
```
