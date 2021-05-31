# A Tour of LazySets

```@meta
CurrentModule = LazySets
```

LazySets is a library for set-based computations in Euclidean space.
The library offers both concrete and lazy set representations,
where the latter stands for the ability to delay set computations until they are
needed. The choice of the programming language Julia and the accompanying documentation
allow researchers to easily translate set-based algorithms
from mathematics to software in a platform-independent way, while achieving runtime
performance that is comparable to statically compiled languages.
Combining lazy operations in high dimensions and explicit computations in low
dimensions, LazySets can be applied to solve complex, large-scale problems.
The library can handle the most common operations between sets, including
[Minkowski sum](https://juliareach.github.io/LazySets.jl/dev/lib/lazy_operations/MinkowskiSum/),
[Cartesian product](https://juliareach.github.io/LazySets.jl/dev/lib/lazy_operations/CartesianProduct/),
[convex hull](https://juliareach.github.io/LazySets.jl/dev/lib/lazy_operations/ConvexHull/), [intersection](https://juliareach.github.io/LazySets.jl/dev/lib/lazy_operations/Intersection/) and
[union](https://juliareach.github.io/LazySets.jl/dev/lib/lazy_operations/UnionSet/).

The puropose of this section is to give a high-level overview of the libary and
several examples to illustrate basic functionality. Further sections delve into
more advanced uses.

## Creating sets

Sets can be created using the constructor for each set type. Here we define the
two-dimensional halfspace

```math
H = \{ (x, y) ∈ \mathbb{R}^2 : x ≤ y\}.
```

```@example tour
using LazySets

H = HalfSpace([1.0, -1.0], 0.0)
```
The first argument is the normal vector to the halfspace and the second argument
the displacement, written in the standard form ``a⋅x ≤ b``. We can define the same
set using *symbolic variables* upon loading the Julia package [Symbolics.jl](https://symbolics.juliasymbolics.org/dev/):

```@example tour
using Symbolics

var = @variables x y

H′ = HalfSpace(x ≤ y, var) # type H\prime[TAB]
```
The resulting sets are equivalent:

```@example tour
# mathematically corresponds to a double inclusion check: H ⊆ H′ && H′ ⊆ H
isequivalent(H, H′)
```
Higher-dimensional sets can be conveniently defined using array notation. Here we
create the halfspace

```math
R = \left\{ x ∈ \mathbb{R}^{10} : \sum_{i=1}^{10} i x_i ≤ 10 \right\}.
```
```@example tour
var = @variables x[1:10]

expr = sum(i*x[i] for i in 1:10) ≤ 10.0
R = HalfSpace(expr, var)
```

This asks LazySets to create two intervals, then take their Cartesian product.

```@example tour
a = Interval(0, 1)
b = Interval(2, 3)

c = a × b
```

In LazySets, intervals are implemented as a thin wrapper around those
available in  [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl).

Use the `convert` method to represent the set in hyperrectangular form (i.e. center and radius vectors):

```@example tour
H = convert(Hyperrectangle, c)
```

Multi-dimensional intervals are of type `IntervalBox`:

```@example tour
using IntervalArithmetic

Hbox = convert(IntervalBox, H)

typeof(Hbox)
```

Hyperrectangles are also special cases of [zonotopes](https://juliareach.github.io/LazySets.jl/dev/lib/sets/Zonotope/).

```@example tour
Z = convert(Zonotope, c)
```

## Operating with sets

The result of **intersecting** two half-spaces is a polyhedron:

```@example tour
var = @variables x y

H1 = HalfSpace(y ≥ x, var)
H2 = HalfSpace(y ≥ -x, var)

P = H1 ∩ H2 # for `∩`, type `\cap[TAB]`
```
Naturally, LazySets can plot the resulting set using the
[Plots.jl](http://docs.juliaplots.org/latest/) library:

```@example tour
using Plots

plot(H1, lab="H1")
plot!(H2, lab="H2")
plot!(P, lab="P = H1 ∩ H2", linewidth=2.0, linestyle=:dash)
```
A polyhedron consists of a finite intersection of half-spaces:

```@example tour
constraints_list(P)
```
The intersection of ``P`` with the triangle ``R`` with vertices ``\{(-1.0, -1.0), (1.0, -1.0), (0.0, 1.0)\}``
can be constructed with:

```@example tour
R = VPolygon([[-1.0, -1.0], [1.0, -1.0], [0.0, 1.0]])

Q = P ∩ R
```
The resulting set is a **lazy (binary) intersection**: by design, all operations are lazy
by default. Concrete operations can be performed calling methods in lower-case:

```@example tour
Q = intersection(P, R)
```
Alternatively, we can use the `concretize` to evaluate a given expression,
transforming it into a suitable concrete set representation.

```@example tour
concretize(P ∩ R)
```
The result is now an `HPolytope` (instead of an `HPolyhedron`), since polytopes
are special cases of *bounded* polyhedra.

```@example tour
isbounded(P)
```

```@example tour
isbounded(Q)
```

```@example tour
plot!(Q, lab="Q = P ∩ T")
plot!(R, lab="R", linestyle=:dashdot, alpha=.2)
```

LazySets can compute various useful functions.

```@example tour
area(Q)
```
By approximating floating-point numbers as rational numbers, we get an exact
value for the area.

```@example tour
area(rationalize(Q))
```

The vertices can be retrieved using `vertices_list`,
```@example tour
vertices_list(Q)
```
or in iterator-fashion with `vertices`:
```@example tour
vertices(Q)
```
The constraint representation of ``Q`` is obtained similarly.

```@example tour
constraints_list(Q)
```
The constraints list iterator is called with `constraints`.
```@example tour
constraints(Q)
```

The **Minkowski sum** of sets can be computed using `⊕` (type `\oplus[TAB]`):

```@example tour
# ball in the 2-norm of given center and radius
B = Ball2([0.8, 1.0], 0.3)

E = Q ⊕ B
```
Approximating and plotting such sets relies on *support function* techniques,
discussed at length in other sections of this manual.

```@example tour
plot(E, lab="E = Q ⊕ B", legend=:topleft)
plot!(Q, lab="Q")
plot!(B, lab="B", ratio=1., xlims=(-0.5, 1.5))
```
**Linear** and **affine** maps are defined combining the `*` and `+` operators.

```@example tour
M = [0 1; 1 0.]
b = [0.5, -0.5]

X = M * E +  b
```
The set ``X``, and in general, the composition of sets and operations, is again
a LazySet (in this case an operation):

```@example tour
isoperation(X)
```
Therefore, all set operations (and plotting) still work in the same way.

```@example tour
plot!(X, lab="X = M*E + b")
```
**Convex hulls** can be computed with `convex_hull` (concrete), `ConvexHull` (lazy, binary),
or `ConvexHullArray` (lazy, n-ary).

```@example tour
C = ConvexHullArray([E, Q, X])

plot!(C, c=:grey, lab="C = CH(E, Q, X)")
```

The **Cartesian product** of sets is computed using the `×` (type `\times[TAB]`).

```@example tour
X = Interval(0, 1) × Ball1([1.0, 2.0, 3.0], 1.0) × Interval(2, 4)
```

**Set unions** are defined with `∪` (type `\cup[TAB]`).

```@example tour
Y = Ball2([1.0, 2.0, 3.0], 1.0) × LineSegment([0.0, 0.0], [1.0, 1.0])

U = X ∪ Y
```

There are other set operations not mentioned in this section. See the remaining
sections of this manual for further examples.

## Exploring the type hierarchy

Every convex set type in the library inherits from the parametric
abstract type `LazySet{N}`, where `N` is a parameter for the numeric
type (typically, double-precision floating point numbers, `Float64`).
This way one can easily switch between, e.g., floating point
(Float64) and exact (Rational) precision with no additional performance penalty:
At runtime, Julia uses multiple dispatch on N and JIT-compiles into type-specific code.

A way to visualize the type hierarchy is using [AbstractTrees](https://github.com/JuliaCollections/AbstractTrees.jl).
The package provides several utilities for working with tree-like data structures.
The functios [`print_tree`](https://juliacollections.github.io/AbstractTrees.jl/stable/api/#AbstractTrees.print_tree)
prints:

```@example tour
using AbstractTrees

AbstractTrees.children(x::Type) = LazySets.subtypes(x, false)

print_tree(LazySet)
```
One of the key features of LazySets (which, admittedly, may be confusing at first)
is that **sets and set operations subtype `LazySet`*. One of the advantages of
these two mathematical concepts to be "just types" is to conveniently compose
representations and operations (existing and user-created ones). An overview
of only the set operations can be obtained like so:

```@example tour
filter(isoperationtype, LazySets.subtypes(LazySet, true))
```
The total amount of available representations is larger (if we had used
`subtypes(LazySet)`, that would only show the direct subtypes of `LazySet`).

```@example tour
length(LazySets.subtypes(LazySet, true))
```

Visualizing the subtypes of a specific set interface can be done similarly.
Consider the class of zonotopic sets, `AbstractZonotope`,
which are those that can be represented as

```math
Z = \left\{ c + ∑_{i=1}^p ξ_i g_i,~~ ξ_i \in [-1, 1]~~ ∀ i = 1,…, p \right\},
```
where ``c \in \mathbb{R}^n`` is called the center and the vectors
``\{g_i\}_{i=1}^p``, ``g_i \in \mathbb{R}^n``, are called the generators.

```@example tour
print_tree(AbstractZonotope)
```

Al the sets belonging to the abstract zonotope interface have in common that they
can be characterized by a center and a generators matrix (equivalently,
as the finite Minkowski sum of line segments, or the the image of a unit infinity-norm
ball in ``\\mathbb{R}^n`` by an affine transformation). Hence, new set types
defining a few interface functions can make use of all the available functionality
that is implemented by generic algorithms, which work at the interface level.

## Random sampling and splitting

Another useful function is to perform **random sampling** from a given set.

```@example tour
# returns a vector of samples (vectors)
S = sample(Q, 300, include_vertices=true)

plot(Q, lab="Q")

# transform into a vector of singletons (for plotting)
plot!(Singleton.(S), lab="S", c=:magenta)
```
Set membership is computed with the `∈` operator (type `\in[TAB]`):

```@example tour
all(s -> s ∈ Q, S)
```

Here is one way to enclose `Q` with smaller, simpler sets (hyperrectangles).

```@example tour
# compute an enclosing hyperrectangle
Qbox = overapproximate(Q, Hyperrectangle)

# split it into 6 pieces horizontally and 10 pieces vertically
B = split(Qbox, [6, 10])

# find all the boxes which intersect Q
idx = findall(Bi -> !isdisjoint(Bi, Q), B)

plot(B[idx])
plot!(Q, lab="Q")
plot!(Singleton.(S), lab="S", c=:magenta)
```

With a bit of extra work we can also split the set into non axis-aligned parallelotopes:

```@example tour
# compute an enclosing zonotope
Z = overapproximate(Q, Zonotope, OctDirections)

# reduce to order 1
P = overapproximate(Z, HParallelotope)

# split along the first (resp. second) generator three (resp. four) times
Zs = split(convert(Zonotope, P), 1:2, [3, 4])

# only keep those zonotopes intersecting Q
idx = findall(Zsi -> !isdisjoint(Zsi, Q), Zs)

plot(Zs[idx])
plot!(Q, lab="Q")
plot!(Singleton.(S), lab="S", c=:magenta)
```


## Polyhedral approximations

In this library we consider representations of closed convex sets in the usual
sense from convex geometry: A set is closed if it contains all its limit points.
A set ``S`` is convex if for any ``m`` points ``v_j ∈ S`` and ``m`` non-negative
numbers ``λ_j`` that sum up to ``1`` we have that ``∑_j λ_j v_j ∈ S`` as well.
Alternatively, a closed convex set is an intersection of (possibly infinitely many)
closed half-spaces.

One of the central functions are `support_function(d, X)` (which admits the alias `ρ(d, X)` (type `\rho[TAB]`))
and `support_vector` (which admits the alias `σ` (type `\sigma[TAB]`)). The support function
of a set ``X`` along direction ``d ∈ \mathbb{R}^n`` corresponds to the (signed) distance
of the maximum element in `X`` along direction ``d``. The result of the support vector
(among possibly infinitely many) that are maximizers of the support function.


## Specialization

One of the key features of LazySets is specialization.
On a more technical level, the library aims at a high level of optimization
to get the best possible performance given the following two restrictions:
the set type and the set dimension.

The submodule `LazySets.Arrays` exports the "one-hot vector" `SingleEntryVector`
that can be used, e.g. to represent polyhedra with axis-aligned half-spaces
as information on the *type*. If `H` is a (non-flat) 100-dimensional hyperrectangle,
it has ``2^{100}`` vertices. Computing the support function of `M*X` for any square
matrix `M` along the canonical direction `e_{50} = [0, 0, …, 1, …, 0]` takes
no more than `10us` in a modern laptop using LazySets.

```@example tour
using LazySets.Arrays: SingleEntryVector

using BenchmarkTools

d = SingleEntryVector(50, 100, 1.0)
M = rand(100, 100)
H = rand(Hyperrectangle, dim=100)
@btime ρ($d, $M * $H)
```

## The Lazy paradigm

For the several set representations defined, such as polyhedra in constraint
and in vertex representation, ellipsoids, balls in different norms, and specific
classes of polyhedra (such as intervals, zonotopes or hyperrectangles),
specialized methods apply. Applications that use LazySets can explore different
approaches with minimal changes in their code. Conversion between set
representations as well as overapproximation or underapproximation algorithms are available.

Set operations can be performed in two possible (complementary) ways:
concretely or lazily, where concrete computation means that the result is a set,
and in the lazy computation the result is an object that represents the operation
between the given sets. For instance, consider the linear map transformation, that is, to compute
``Y = \{y : y = Ax \textrm{ for some } x ∈ X\}``. In LazySets,
`linear_map(A, X)` returns a *concrete* set representation of ``Y``.
The algorithm that is actually used, as well as the type of ``Y``, depends on the types of
its arguments and the dimension of X (this is the *multiple dispatch* paradigm).

For example, if `X` is a polygon in vertex representation (`VPolygon`)
then `linear_map(M, X)` applies `M` to each vertex of `X`. However, if `X` is a 30-dimensional
polyhedron in halfspace representation (`HPolyhedron`), then the halfspace representation
is used if `M` is invertible -- if not, and if `MX` doesn't fall into any special case,
the vertex representation has to be computed using the polyhedra manipualation library
`Polyhedra.jl` (the cost of computing the vertex representation increases
exponentially with the dimension, so computing with concrete polyhedra is usually avoided
in high dimensions). Observe that if we are interested in `P(MX)` for some projection matrix `P`,
that operation *can* be performed efficiently using support functions, i.e. without actually computing `MX`.
On the other hand, there are special cases of polyhedra for which concrete operations
can be performed efficiently. For instance, if `X` is a hyperrectangle (or a zonotope),
the resulting set is a zonotope and the computation is done efficiently.

On the other hand, `LinearMap(A, B)`, or simply `A * X`, computes the lazy linear map,
which is just a new object which "holds" the linear map computation until it is
actually needed. In other words, LinearMap(A, X) can be used to reason about the
linear map even if X is very high-dimensional, since this command just builds an
object representing the linear map of `X` by `A`.
In Julia, Unicode symbols such as `A * X`, `X ⊕ Y`, `X ⊖ Y`, `X × Y`, all default
to lazy operations by design. Unicode is input using the latex name of the symbol followed by the TAB key,
such as `\oplus[TAB]` for the (lazy) Minkowski sum ⊕.

Finally, one can combine lazy set operations to build lazy expressions that
represent several operations between sets, such as eg. `Q = (Z ⊕ A*X) × T`.
By means of the basic tools of convex geometry, useful information about `Q` can be obtained without actually
computing the linear map, minkowski sum and cartesian product in the above computation.
If the computation only involves querying information about `Q` in a restricted
number of directions, as it is often the case in applications, the lazy approach
can be realized very efficiently!

## How to contribute

Development happens on [github](https://github.com/JuliaReach/LazySets.jl/).
New contributors should follow the links provided in the
[About](https://juliareach.github.io/LazySets.jl/dev/about/) section of this documentation.

## How to cite

When citing `LazySets.jl`, please cite [JuliaReach: a toolbox for set-based reachability
](https://dl.acm.org/doi/10.1145/3302504.3311804).

```
@inproceedings{bogomolov2019juliareach,
  title={JuliaReach: a toolbox for set-based reachability},
  author={Bogomolov, Sergiy and Forets, Marcelo and Frehse, Goran and Potomkin, Kostiantyn and Schilling, Christian},
  booktitle={Proceedings of the 22nd ACM International Conference on Hybrid Systems: Computation and Control},
  pages={39--44},
  year={2019}
}
```
