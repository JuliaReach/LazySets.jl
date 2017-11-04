# LazySets.jl

`LazySets` is a [Julia](http://julialang.org) package for calculus with
convex sets.

The aim is to provide a scalable library for solving complex set-based
problems, such as those encountered in [differential inclusions](https://en.wikipedia.org/wiki/Differential_inclusion)
or reachability analysis techniques in the domain of [formal verification](https://en.wikipedia.org/wiki/Formal_verification).
Typically, one is confronted with a set-based recurrence with a given initial set and/or
input sets, and for visualization purposes the final result has to be obtained through an
adequate projection onto low-dimensions. This library implements types to construct
the set formulas with the usual mathematical notation, and methods to efficiently
and accurately approximate the projection in low-dimensions.

```@contents
Pages = ["index.md"]
```

## Introduction

The strategy consists of defining lazy (i.e. symbolic) representations that are used
to write set-based formulas. This provides an exact but abstract formulation and may involve any common convex set
class or operation between sets. Then, concrete information is obtained through querying specific directions.
More precisely, each concrete subtype $\mathcal{X}$ of the abstract type `LazySet`,
exports a method to calculate its support vector $\sigma(d, \mathcal{X})$ in a given (arbitrary) direction
$d \in \mathbb{R}^n$. Representing sets exactly but lazily has the advantage of being
able to perform only the required operations on-demand.

For very long computations (e.g. set-based recurrences with tens of thousands of elements),
it is useful to combine both lazy and concrete representations such as polyhedral approximations.
All this is easy to do with `LazySets`. Moreover, there is a specialized module for handling
arrays of two-dimensional projections using Cartesian decomposition techniques. The projection
can be taken to the desired precision using an iterative refinement method.

## Example

Let $\mathcal{X}_0 \subset \mathbb{R}^{1000}$ be the Euclidean ball of center
$(1, \ldots, 1)$ and radius $0.1$ in dimension $n=1000$. Suppose that given a real matrix
$A \in \mathbb{R}^{1000 \times 1000}$ we are interested in the equation:

$\mathcal{Y} = CH(e^{A δ} \mathcal{X}_0 ⊕ δ B\mathcal{U}, \mathcal{X}_0),$

where $CH$ is the convex hull operator, $⊕$ denotes Minkowski sum,
$\mathcal{U}$ is a ball in the infinity norm centered at zero and radius $1.2$,
and $B$ is a linear map of the appropriate dimensions.

For concreteness, let's take $A$ to be a random matrix with probability $1\%$ of any entry being nonzero.
Let's suppose that the input set $\mathcal{U}$ is two-dimensional, and that the linear map $B$ is random.
Using `LazySets` we can define this problem as follows:

```julia
using LazySets
A = sprandn(1000, 1000, 0.01)
δ = 0.1
X0 = Ball2(ones(1000), 0.1)
B = randn(1000, 2)
U = BallInf(zeros(2), 1.2)
```

The `@time` macro reveals that building $\mathcal{Y}$ with
`LazySets` is instantaneous:

```julia
@time Y = CH(SparseMatrixExp(A * δ) * X0 + δ * B * U, X0);
0.000022 seconds (13 allocations: 16.094 KiB)
```

By asking the concrete type of `Y`, we see that this object is a convex hull type,
parameterized by the types of its arguments, corresponding to the mathematical formulation:

```julia
julia> typeof(Y)
LazySets.ConvexHull{LazySets.MinkowskiSum{LazySets.ExponentialMap{LazySets.Ball2},
LazySets.LinearMap{LazySets.BallInf}},LazySets.Ball2}
```

Now suppose that we are interested in observing the projection of $\mathcal{Y}$ onto
the variables (1, 500). First we define the $2×1000$ projection matrix and apply it to $\mathcal{Y}$
as a linear map (i.e. from the left). Second, we use the `overapproximate` method:

```julia
proj_mat = [[1. zeros(1, 999)]; [zeros(1, 499) 1. zeros(1, 500)]]
@time res = Approximations.overapproximate(proj_mat * Y);
0.064034 seconds (1.12 k allocations: 7.691 MiB)
```

We have calculated a box overapproximation of the exact projection onto the (1, 500) plane.
Notice that it takes about 0.064 seconds for the whole operation, allocating less than
10MB or RAM. Let us note that if the set operations where done explicitly,
this would be much (!) slower. For instance, already the explicit computation of the matrix exponential would
have costed 10x more, and allocated around 300MB. For even higher $n$,
you'll probably run out of RAM! But this is doable with `LazySets` because the *action* of the matrix exponential over the set
is being computed, evaluated only along the directions of interest.
Similar comments apply to the Minkowski sums above.

We can visualize the result using `plot`, as shown below (left-most plot).

![assets/example_ch.png](assets/example_ch.png)

In the second and third plots, we have used a refined method that allows to specify a prescribed accuracy
for the projection (in terms of [Hausdorff distance](https://en.wikipedia.org/wiki/Hausdorff_distance)).
It can be passed as a second argument to `overapproximate`. 

|Error tol.|time (s)|memory (MB)|
|------|------|------|
|∞ (no refinement)|0.022|5.27|
|1e-1|0.051|7.91|
|1e-3|0.17|30.3|

This table shows the runtime and memory consumption for different error tolerances,
and the results are shown in three plots of above, from left to right. When passing
to a smaller tolerance, the corners connecting edges are more "rounded", and naturally,
for increasing accuracy it becomes more expensive computationally, since more
support vectors have to be evaluated.

## Features

The core functionality of `LazySets` is:

- Lazy (i.e. symbolic) types for several classes of convex sets such as
  balls in different norms, polygons in constraint or vertex representation,
  special types such as lines and linear constraints, hyperrectangles, and
  high-dimensional polyhedra.
- Most commonly used set operations, e.g. Minkowski sum, Cartesian product,
  convex hull and interval hull approximations. Moreover, lazy linear maps and
  lazy exponential maps are also provided.

On top of the previous basic type representations and operations, `LazySets`
can be used to:

- Efficiently evaluate the support vector of nested lazy sets using parametrized `LazySet` arrays.
- Cartesian decomposition of lazy sets using two-dimensional projections.
- Fast overapproximation of an exact set using a polyhedral
  approximation, to the desired accuracy.
- Extensive visualization capabilities through Julia's
  [Plots.jl](http://docs.juliaplots.org/latest/) framework.

## Manual Outline

```@contents
Pages = [
    "man/getting_started.md",
    "man/support_vectors.md",
    "man/polyhedral_approximations.md",
    "man/fast_2d_LPs.md",
    "man/iterative_refinement.md"
]
Depth = 2
```

## Library Outline

```@contents
Pages = [
    "lib/representations.md",
    "lib/operations.md",
    "lib/approximations.md",
    "lib/utils.md"
]
Depth = 2
```
