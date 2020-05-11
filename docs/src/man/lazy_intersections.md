# Lazy Intersections

In this section we illustrate the use of lazy intersection in `LazySets`. We
will use the ellipsoid set type.

```@contents
Pages = ["lazy_intersections.md"]
Depth = 3
```

An ellipsoid ``E`` can be created by giving its center ``c`` and its shape matrix ``Q``,
which should be positive definite, i.e. its eigenvalues must be positive.
Mathematically, it is the set

```math
    E = \{ x ∈ \mathbb{R}^n : (x-c)Q^{-1}(x-c) ≤ 1\}.
```

Let's make two rotated ellipsoids and plot them in the same pair of axes.

```@example ellipsoids
using Plots, LazySets

E₁ = Ellipsoid(zeros(2), [1 0; 0 2.])
E₂ = Ellipsoid(ones(2), [2 0; 0 1.])

pell = plot(E₁, aspectratio=1)
pell = plot!(pell, E₂)
```

!!! note
    The accuracy to which this set is plotted can be controlled by passing a numerical argument as in
    `plot(E₁, 1e-3, aspectratio=1)`. Here, `1e-3` stands for an upper-bound on the error, measured in terms of
    the [Hausdorff distance](https://en.wikipedia.org/wiki/Hausdorff_distance) between the ellipsoid and
    the polygonal overapproximation which is actually computed for display.

Now let's take the *lazy* intersection of the ellipses:


```@example ellipsoids
Z = E₁ ∩ E₂
typeof(Z)
```

On the other hand, the *concrete* intersection of sets, called `intersection` in
`LazySets`, is not yet available for ellipsoids:

```@example ellipsoids
hasmethod(intersection, Tuple{typeof(E₁), typeof(E₂)})
```

So how can we work with the intersection of the ellipsoids?

One way is to overapproximate them by polygons (or polytopes in higher dims) and
then take their intersection, because this function is defined, whose return
type is again a `HPolytope`:

```@example ellipsoids
hasmethod(intersection, Tuple{HPolytope{Float64}, HPolytope{Float64}})
```

```@example ellipsoids
import LazySets.Approximations.overapproximate

# the parameter epsilon controls the accuracy of the iterative refinement,
# with respect to the Hausdorff distance
H₁(ε) = overapproximate(E₁, HPolygon, ε)
H₂(ε) = overapproximate(E₂, HPolygon, ε)

# using the concrete hpolytope-hpolytope intersection here
Hint(ε) = intersection(convert.(HPolytope, [H₁(ε), H₂(ε)])...);
```

```@example ellipsoids
pell = plot(E₁, aspectratio=1)
pell = plot!(pell, E₂)
pεsmaller = plot!(pell, convert(HPolygon, Hint(0.5)), alpha=.4)

pell = plot(E₁, aspectratio=1)
pell = plot!(pell, E₂)
pεbigger = plot!(pell, convert(HPolygon, Hint(0.05)), alpha=.4)

plot(pεsmaller, pεbigger, layout=(1, 2))
```

Note how dividing the $\varepsilon$ threshold by 10 makes the polygonal
overapproximation of the intersection tighter.

Yet another approach is to directly query the directions of the *lazy* intersection
`E₁ ∩ E₂`. We can overapproximate using template directions, such as a box,
an octagon, or other.

The idea behind the template overapproximation method is to use the property that
the support function of the intersection of two convex sets is upper bounded by
the min of the support function of each set. We can see in the following experiments that
the resulting set is quite tight.

```@example ellipsoids
using Polyhedra

# overapproximate the lazy intersection using a box
Xbox = overapproximate(E₁ ∩ E₂, BoxDirections(2))

# overapproximate the lazy intersection using octagonal directions
Xoct = overapproximate(E₁ ∩ E₂, OctDirections(2))

pell = plot(E₁, aspectratio=1)
pell = plot!(pell, E₂)
pbox = plot!(pell, Xbox, alpha=.4)

pell = plot(E₁, aspectratio=1)
pell = plot!(pell, E₂)
poct = plot!(pell, Xoct, alpha=.4)

plot(pbox, poct, layout=(1, 2))
```

Using support function evaluations over a set of fixed directions is in general
more efficient than iterative refinement, but the drawback is that one does not have control on the
overapproximation error. Moreover, iterative refinement is currently only available in two dimensions,
but overapproximation with template directions can be used in any dimension.

Let's time it!

```@example ellipsoids
using BenchmarkTools

@btime overapproximate($E₁ ∩ $E₂, BoxDirections(2))
@btime overapproximate($E₁ ∩ $E₂, OctDirections(2));
```

We can work with higher dimensional ellipsoids as well:

```@example ellipsoids
using LinearAlgebra

# a random ellipsoid in n-dimensions
function rand_ellipsoid(n)
    A = rand(n,n)
    Q = (A+transpose(A))/2 + n * I
    Ellipsoid(rand(n), Q)
end;
```

```@example ellipsoids
for n in [2, 5, 50, 100]
    println("\nn = $n\n")
    global E₁, E₂ = rand_ellipsoid(n), rand_ellipsoid(n)

    # overapproximate the lazy intersection using an n-dimensional box
    @btime overapproximate($E₁ ∩ $E₂, BoxDirections($n))
    
    # overapproximate the lazy intersection using octagonal directions in R^n
    @btime overapproximate($E₁ ∩ $E₂, OctDirections($n))
end;
```
