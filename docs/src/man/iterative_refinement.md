# Iterative Refinement

This section of the manual describes an approximation method for an arbitrary
two-dimensional convex set ``S`` and a given error bound ``ɛ`` using support
vectors.
The basic idea is to add new supporting directions whenever the approximation
error is still bigger than ``ɛ``.

```@contents
Pages = ["iterative_refinement.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets.Approximations
DocTestSetup = quote
    using LazySets, Plots, LazySets.Approximations
end
```

## Local approximations

The polygonal approximation of an arbitrary lazy convex set `S` is represented
by a list of local approximations or refinements.
More precisely, a *local approximation* is a triple ``(p_1, p_2, q)``, where:

- ``p_1`` and ``p_2`` belong to ``S``
- the segments ``(p_1 q)`` and ``(p_2 q)`` belong to support lines of ``S``

Since ``S`` is assumed to be convex, the segment ``(p_1 p_2)`` is inside ``S``.
Taking each support line ``(p_1 q)`` of a given list of local approximations of
``S``, we can build a polygon in constraint representation that overapproximates
`S`.

The type `LocalApproximation{N}` implements a local approximation; it is
parametric in the numeric type `N`, and also contains additional information
regarding the quality of the approximation:
The `refinable` field is a boolean that is `true` whenever the approximation can
be improved, and `err` is an upper bound on the exact Hausdorff distance of the
approximation with respect to the exact set `S`.

Given the unit ball in the 2-norm, below we plot the local approximation along
the East and North directions.

```@example example_iterative_refinement
using LazySets, Plots, LazySets.Approximations

b = Ball2(zeros(2), 1.)

plot(b, 1e-3, aspectratio=1, alpha=0.3)

plot!(Singleton([1.0, 0.0]), annotations=(1.1, 0.1, text("p1")), color="green")
plot!(Singleton([0.0, 1.0]), annotations=(0.1, 1.1, text("p2")), color="green")
plot!(Singleton([1.0, 1.0]), annotations=(1.09, 1.1, text("q")))
plot!(Singleton([0.0, 0.0]), annotations=(0.1, 0.0, text("0")), color="green")
plot!(annotations=(1.4, 0.1, text("d1")))
plot!(annotations=(0.1, 1.4, text("d2")))
plot!(annotations=(0.75, 0.8, text("ndir")))

plot!(x->x, x->1., -0.8, 1.3, line=1, color="black", linestyle=:dash)
plot!(x->1., x->x, -0.8, 1.3, line=1, color="black", linestyle=:dash)
plot!(x->x+1, x->0., 0.0, 0.4, line=1, color="red", linestyle=:solid, arrow=true)
plot!(x->0., x->x+1, 0.0, 0.4, line=1, color="red", linestyle=:solid, arrow=true)
plot!(x->-x, x->x+1, -1.2, .2, line=1., color="black", linestyle=:dashdot)
plot!(x->x+.6, x->x+.6, -.1, .08, line=1, color="red", linestyle=:solid, arrow=true)
```

We can instantiate and append this approximation to a fresh
`PolygonalOverapproximation` object, which is a type that wraps a set and a list
of `LocalApproximation`s.
The approximation is refinable, since it can be "split" along `ndir`, where
`ndir` is the direction normal to the line ``(p_1 p_2)`` (shown dash-dotted in
the figure), providing two approximations which are closer to the given set in
Hausdorff distance.


```@example example_iterative_refinement
import LazySets.Approximations:PolygonalOverapproximation, addapproximation!

Ω = PolygonalOverapproximation{Float64}(b)
p1, d1, p2, d2 = [1.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 1.0]
approx_EAST_NORTH = addapproximation!(Ω, p1, d1, p2, d2)

approx_EAST_NORTH.refinable
```

The associated error is ``\sqrt{2}-1≈0.414213``, which is the distance between
the point ``q`` and the intersection between the line ``(0 q)`` and the circle.
Actually this point corresponds to the support vector of the set `b` along
`ndir`.

```@example example_iterative_refinement
approx_EAST_NORTH.err
```

The refined approximation is computed next.

## Refinement

Basically, the refinement step consists of splitting the local approximation
``(p_1, p_2, q)`` into two local approximations ``(p_1, s, q')`` and
``(s, p_2, q'')``, where `s` is the support vector of ``S`` along `ndir`.

To illustrate this, first let's add the remaining three approximations to `Ω`
along the canonical directions, to build a box overapproximation of `b`.

```@example example_iterative_refinement
import LazySets.Approximations: refine, tohrep

plot(b, 1e-3, aspectratio=1, alpha=0.3)

# initialize box directions
DIR_EAST, DIR_NORTH, DIR_WEST, DIR_SOUTH = [1., 0.], [0., 1.], [-1., 0.], [0., -1.]
pe, pn, pw, ps = σ(DIR_EAST, b), σ(DIR_NORTH, b), σ(DIR_WEST, b), σ(DIR_SOUTH, b)

addapproximation!(Ω, pn, DIR_NORTH, pw, DIR_WEST)
addapproximation!(Ω, pw, DIR_WEST, ps, DIR_SOUTH)
addapproximation!(Ω, ps, DIR_SOUTH, pe, DIR_EAST)

plot!(tohrep(Ω), alpha=0.2, color="orange")
```

Next we refine the first approximation of the list.

```@example example_iterative_refinement
(r1, r2) = refine(Ω, 1)
Ω.approx_list[1] = r1
insert!(Ω.approx_list, 2, r2)

plot(b, 1e-3, aspectratio=1, alpha=0.3)
plot!(tohrep(Ω), alpha=0.2, color="orange")
```

We call `r1` and `r2` the right and left approximations respectively, since they
are saved in counter-clockwise order.
We can check that the first two approximations are still refinable.

```@example example_iterative_refinement
Ω.approx_list[1].refinable,  Ω.approx_list[2].refinable
```

Hence, we can make again a refinement of that approximation.

```@example example_iterative_refinement
(r1, r2) = refine(Ω, 1)
Ω.approx_list[1] = r1
insert!(Ω.approx_list, 2, r2)

plot(b, 1e-3, aspectratio=1, alpha=0.3)
plot!(tohrep(Ω), alpha=0.2, color="orange")
```

The criterion for an approximation being refinable is that we can properly
define a normal direction `ndir`.
This boils down to checking for the following "degenerate" cases:

1. ``p_1`` and ``p_2`` overlap.
2. ``p_1`` and ``q`` overlap.
1. ``p_2`` and ``q`` overlap.

Moreover, we include the condition `approx_error > TOL` where `TOL` is the floating
point epsilon in the given numerical precision.

## Algorithm

Having presented the individual steps, we give the pseudocode of the iterative
refinement algorithm, see `approximate(S, ε)`.

The algorithm consists of the following steps:

1. *Initialization*. The approximation is initialized with box directions,
   i.e. it starts with four `LocalApproximation` objects. Let `i=1`.
2. *Refinement loop*. If the local approximation at index `i` has an error
   greater than the threshold `ε`, then refine. Otherwise, increment `i <- i+1`.
3. *Redundancy check*. Insert the refined right approximation at position `i`,
   and check whether the left approximation is redundant or not with respect to
   the one at position `i+1`. Checking for redundancy amounts to checking for
   overlap of both `p1` and `q`. Then, either substitute at `i+1` or insert
   (keeping the approximation at `i+1`) depending on the redundancy check.
4. *Stopping criterion*. Terminate if the index `i` exceeds the current length
   of the approximations list; otherwise continue with step 2.

Observe that the algorithm finishes when all approximations are such that
their associated error is smaller than `ε`, hence the Hausdorff distance between
`S` and its polygonal overapproximation is no greater than `ε`.

## Example

As a final example consider the iterative refinement of the ball `b` for
different values of the approximation threshold `ε`.

```@example example_iterative_refinement
import LazySets.Approximations:overapproximate, approximate

p0 = plot(b, 1e-6, aspectratio=1)
p1 = plot!(p0, overapproximate(b, 1.), alpha=0.4, aspectratio=1)

p0 = plot(b, 1e-6, aspectratio=1)
p2 = plot!(p0, overapproximate(b, 0.1), alpha=0.4, aspectratio=1)

p0 = plot(b, 1e-6, aspectratio=1)
p3 = plot!(p0, overapproximate(b, 0.01), alpha=0.4, aspectratio=1)

plot(p1, p2, p3, layout=(1, 3))
```

We can check that the error is getting smaller with `ε` in each case:

```@example example_iterative_refinement
f = x -> (minimum(x), maximum(x))
g = ε -> f([ai.err for ai in approximate(b, ε).approx_list])
g(1.), g(0.1), g(0.01)
```

Meanwhile, the number of constraints of the polygonal overapproximation
increases, in this example by a power of 2 when the error is divided by a factor
10.

```@example example_iterative_refinement
h = ε ->  length(approximate(b, ε).approx_list)
h(1.), h(0.1), h(0.01)
```

!!! note
    Actually, the plotting function for an arbitrary `LazySet` `plot(...)`,
    called *recipe* in the context of
    [Plots.jl](https://github.com/JuliaPlots/Plots.jl), is such that it receives
    a numeric argument `ε` and the routine itself calls `overapproximate`.
    However, some sets such as abstract polygons have their own plotting recipe
    hence do not require the error threshold, since they are plotted exactly as
    the convex hull of their vertices.
