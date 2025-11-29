# A Tour of LazySets

```@meta
CurrentModule = LazySets
```

LazySets is a library for set-based computations in Euclidean space.
The library offers both concrete and lazy set representations,
where the latter stands for the ability to delay set computations until they are
needed. The choice of the programming language [Julia](https://julialang.org/)
and the accompanying documentation
allow researchers to easily translate set-based algorithms
from mathematics to software in a platform-independent way, while achieving
runtime performance that is comparable to statically compiled languages.
Combining lazy operations in high dimensions and explicit computations in low
dimensions, LazySets can be applied to solve complex, large-scale problems.
The library can handle the most common operations between sets, including
[Minkowski sum](https://juliareach.github.io/LazySets.jl/dev/lib/lazy_operations/MinkowskiSum/),
[Cartesian product](https://juliareach.github.io/LazySets.jl/dev/lib/lazy_operations/CartesianProduct/),
[convex hull](https://juliareach.github.io/LazySets.jl/dev/lib/lazy_operations/ConvexHull/),
[intersection](https://juliareach.github.io/LazySets.jl/dev/lib/lazy_operations/Intersection/)
and [union](https://juliareach.github.io/LazySets.jl/dev/lib/lazy_operations/UnionSet/).

The purpose of this section is to give a high-level overview of the library and
several examples to illustrate basic functionality. Further sections delve into
more advanced uses.


## Creating sets

Sets can be created using the constructor for each set type. Here we define the
two-dimensional half-space

```math
H = \{ (x, y) ∈ ℝ^2 : x ≤ y\} = \{ (x, y) ∈ ℝ^2 : x - y ≤ 0\}.
```

```@example tour
using LazySets

H = HalfSpace([1.0, -1.0], 0.0)
```

The first argument is the normal vector to the half-space and the second
argument the displacement, written in the standard form ``a⋅x ≤ b``. We can
define the same set using *symbolic variables* with the help of the Julia
package [Symbolics.jl](https://symbolics.juliasymbolics.org/dev/):

```@example tour
using Symbolics

var = @variables x y

H′ = HalfSpace(x ≤ y, var)  # for `′`, type H\prime<tab>
```

We can use LazySets to check that `H` and `H′` represent the same sets
(mathematically the set equivalence ``H = H'`` corresponds to a double inclusion
check ``H ⊆ H' ∧ H' ⊆ H``):

```@example tour
isequivalent(H, H′)
```

High-dimensional sets can be conveniently defined using array notation. Here we
create the nine-dimensional half-space

```math
H9 = \left\{ x ∈ ℝ^{9} : ∑_{i=1}^{9} i x_i ≤ 10 \right\}.
```

```@example tour
var = @variables x[1:9]

expr = sum(i*x[i] for i in 1:9) ≤ 10.0
H9 = HalfSpace(expr, var)
```

Next we ask LazySets to create two intervals ``[0, 1]`` and ``[2, 3]`` and then
take their Cartesian product.

```@example tour
A = Interval(0, 1)
B = Interval(2, 3)

C = A × B
```

In LazySets, intervals are implemented as thin wrappers around those available
in [IntervalArithmetic.jl](https://juliaintervals.github.io/pages/packages/intervalarithmetic/).
Use the `convert` function to represent the set as a `Hyperrectangle`, which is
a dedicated set type that represents a hyperrectangle with a center and a radius
vector:

```@example tour
H = convert(Hyperrectangle, C)
```

IntervalArithmetic.jl represents multi-dimensional intervals with the type
`IntervalBox`, to which we can also `convert`:

```@example tour
using IntervalArithmetic: IntervalBox

Hbox = convert(IntervalBox, H)

typeof(Hbox)
```

Hyperrectangles are a special subclass of
[zonotopes](https://juliareach.github.io/LazySets.jl/dev/lib/sets/Zonotope/).

```@example tour
Z = convert(Zonotope, C)
```

There are many other set representations available in this library. They are
introduced in further sections below.


## Operating with sets

The result of **intersecting** two or more half-spaces is a polyhedron:

```@example tour
var = @variables x y

H1 = HalfSpace(y ≥ x, var)
H2 = HalfSpace(y ≥ -x, var)

P = H1 ∩ H2  # for `∩`, type `\cap<tab>`
```

Naturally, LazySets can plot the resulting set using the
[Plots.jl](http://docs.juliaplots.org/latest/) library:

```@example tour
using Plots

plot(H1, lab="H1")
plot!(H2, lab="H2")
plot!(P, lab="P = H1 ∩ H2", linewidth=2.0, linestyle=:dash)
```

A polyhedron consists of a finite intersection of half-spaces, which are also
called linear constraints:

```@example tour
constraints_list(P)
```

The intersection of ``P`` with the triangle ``R`` with vertices

``\{(-1.0, -1.0), (1.0, -1.0), (0.0, 1.0)\}`` can be constructed with:

```@example tour
R = VPolygon([[-1.0, -1.0], [1.0, -1.0], [0.0, 1.0]])

Q = P ∩ R
```

```@example tour
plot!(Q, lab="Q = P ∩ R")
plot!(R, lab="R", linestyle=:dashdot, alpha=.2)
```

The resulting set `Q` is a **lazy (binary) intersection**: by design, all
operations in LazySets are lazy by default. Concrete operations can be performed
by calling corresponding lower-case functions:

```@example tour
Q = intersection(P, R)
```

Alternatively, we can use the `concretize` function to transform a lazy
set representation into a suitable concrete set representation.

```@example tour
concretize(P ∩ R)
```

The result is an `HPolytope` (instead of an `HPolyhedron`), since polytopes are
special cases of *bounded* polyhedra.

```@example tour
isbounded(P)
```

```@example tour
isbounded(Q)
```

LazySets can compute various useful operations. For example, it can compute the
area of a two-dimensional bounded set:

```@example tour
area(Q)
```

When using rational numbers instead of floating-point numbers, we get an exact
value for the area:

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

The constraint iterator is called with `constraints`.

```@example tour
constraints(Q)
```

The **Minkowski sum** of two sets can be computed using `⊕` (type `\oplus<tab>`)
(or alternatively just `+` for convenience):

```@example tour
# ball in the 2-norm of given center and radius
B = Ball2([0.8, 1.0], 0.3)

E = Q ⊕ B
```

Approximating and plotting such sets relies on techniques based on the
*support function*, discussed at length in other sections of this manual.

```@example tour
plot(E, lab="E = Q ⊕ B", legend=:bottomright)
plot!(Q, lab="Q")
plot!(B, lab="B", ratio=1., xlims=(-0.5, 1.5))
```

**Linear** and **affine** maps are defined combining the `*` and `+` operators.

```@example tour
M = [0 1; 1 0.]
b = [0.5, -0.5]

X = M * E + b
```

The set ``X``, and in general, the composition of sets and operations, is again
a set (in this case an operation).

```@example tour
isoperation(X)
```

Therefore, all set operations (and plotting) still work in the same way.

```@example tour
plot!(X, lab="X = M*E + b", legend=:topright)
```

**Convex hulls** can be computed with `convex_hull` (concrete), `ConvexHull`
(lazy, binary), or `ConvexHullArray` (lazy, n-ary).

```@example tour
C = ConvexHullArray([E, Q, X])

plot!(C, c=:grey, lab="C = CH(E, Q, X)")
```

The **Cartesian product** of sets is defined with `×` (type `\times<tab>`) (or
alternatively just `*` for convenience).

```@example tour
X = Interval(0, 1) × Ball1([1.0, 2.0, 3.0], 1.0) × Interval(2, 4)
```

**Set unions** are defined with `∪` (type `\cup<tab>`).

```@example tour
Y = Ball2([1.0, 2.0, 3.0], 1.0) × LineSegment([0.0, 0.0], [1.0, 1.0])

U = X ∪ Y
```

There are other set operations not mentioned in this section. See the remaining
sections of this manual for further examples.


## Exploring the type hierarchy

Every set type in the library inherits from the parametric abstract type
`LazySet{N}`, where `N` is a parameter for the numeric type (typically,
double-precision floating point numbers, `Float64`).
This way one can easily choose between, e.g., floating point (`Float64`) and
exact (`Rational`) precision with no additional performance penalty:
At runtime, Julia uses multiple dispatch on `N` and JIT-compiles into
type-specific code.

Since the `LazySets` type hierarchy is rather involved, a way to visualize it is
to use [AbstractTrees.jl](https://juliacollections.github.io/AbstractTrees.jl/stable/).
The package provides several utilities for working with tree-like data structures.
The function
[`print_tree`](https://juliacollections.github.io/AbstractTrees.jl/stable/api/#AbstractTrees.print_tree)
provides a very detailed answer.

```@example tour
using AbstractTrees

AbstractTrees.children(x::Type) = LazySets.subtypes(x, false)

print_tree(LazySet)
```

Since the list does not fit into the default size, some types are hidden.
LazySets has several set interfaces whose names start with `Abstract`.
Visualizing the subtypes of a specific set interface can be done similarly.
Consider the class of zonotopic sets, `AbstractZonotope`, which are those that
can be represented as

```math
Z = \left\{ x ∈ ℝ^n : x = c + ∑_{i=1}^p ξ_i g_i,~~ ξ_i ∈ [-1, 1]~~ ∀ i = 1,…, p \right\},
```

where ``c ∈ ℝ^n`` is called the center and the vectors
``\{g_i\}_{i=1}^p``, ``g_i ∈ ℝ^n``, are called the generators.

```@example tour
print_tree(AbstractZonotope)
```

All the sets belonging to the abstract zonotope interface have in common that
they can be characterized by a center and a generator matrix (equivalently, as
the finite Minkowski sum of line segments, or as the image of a unit
infinity-norm ball in ``ℝ^n`` by an affine transformation). Hence, new
set types implementing a few interface functions can make use of all the
available functionality that are implemented by generic algorithms working at
the interface level.

Last but not least, we remark that one of the key design choices of LazySets
(which, admittedly, may be confusing at first) is that
**both set representations and set operations subtype `LazySet`**.
In other words, an "operation between sets" and "a set" are on the same footing
in terms of belonging to the same type hierarchy. One of the advantages of these
two mathematical concepts to be "just types" is to conveniently
compose representations and operations (existing and user-created ones).

An overview of only the set operations can be obtained like so:

```@example tour
filter(isoperationtype, LazySets.subtypes(LazySet, true))
```

The total amount of available representations is larger (if we had used
`subtypes(LazySet)`, that would only show the direct subtypes of `LazySet`).

```@example tour
length(LazySets.subtypes(LazySet, true))
```


## Random sampling and splitting

Another useful operation is to perform **random sampling** from a given set.

```@example tour
# returns a vector of samples (vectors)
S = sample(Q, 300, include_vertices=true)

plot(Q, lab="Q")

# transform into a vector of singletons (for plotting)
plot!(Singleton.(S), lab="S", c=:magenta)
```

Set membership is computed with the `∈` operator (type `\in<tab>`):

```@example tour
all(s -> s ∈ Q, S)
```

Below we show a straightforward way to enclose `Q` with smaller, simpler sets
(hyperrectangles in this case).

```@example tour
# compute an enclosing hyperrectangle
Qbox = overapproximate(Q, Hyperrectangle)

# split it into 6 pieces horizontally and 10 pieces vertically
Bs = split(Qbox, [6, 10])

# find all the boxes that intersect Q
idx = findall(Bi -> !isdisjoint(Bi, Q), Bs)

plot(Bs[idx])
plot!(Q, lab="Q")
plot!(Singleton.(S), lab="S", c=:magenta)
```

With a bit of extra work we can also split the set into non axis-aligned
parallelotopes:

```@example tour
# compute an enclosing zonotope
Z = overapproximate(Q, Zonotope, OctDirections)

# reduce to order 1
P = overapproximate(Z, HParallelotope)

# split along the first (resp. second) generator three (resp. four) times
Zs = split(convert(Zonotope, P), [1, 2], [3, 4])

# find all the zonotopes that intersect Q
idx = findall(Zi -> !isdisjoint(Zi, Q), Zs)

plot(Zs[idx])
plot!(Q, lab="Q")
plot!(Singleton.(S), lab="S", c=:magenta)
```


## Polyhedral approximations

In this library we mainly consider representations of closed convex sets in the
usual sense from convex geometry: A set is closed if it contains all its limit
points. A set ``S`` is convex if for any ``m`` points ``v_j ∈ S`` and ``m``
non-negative numbers ``λ_j`` that sum up to ``1`` we have that
``∑_{j=1}^m λ_j v_j ∈ S`` as well. Alternatively, a closed convex set is an
intersection of (possibly infinitely many) closed half-spaces.

One of the central functions are `ρ(d, X)` (type `\rho<tab>`; with the alias
`support_function`) and `σ` (type `\sigma<tab>`; with the alias `support_vector`).
The support function of a set ``X`` along direction ``d ∈ ℝ^n`` corresponds to
the (signed) distance of the maximum element in ``X`` along direction ``d``.
The support vector is a (among possibly infinitely many) maximizer of the support
function.

The support function also characterizes the half-space that tightly covers a set
in a given direction. Hence it can be used to obtain a polyhedral
(over-)approximation of any set. This topic is further covered in other sections
of the manual.


## Specialization

One of the key features of LazySets is specialization.
On a more technical level, the library aims at a high level of optimization
to get the best possible performance given the following two restrictions:
the set type and the set dimension.

The submodule `LazySets.Arrays` exports the "one-hot vector" `SingleEntryVector`
that can be used, e.g., to represent polyhedra with axis-aligned half-spaces
as information on the *type*. If `H` is a (non-flat) 100-dimensional
hyperrectangle, it has ``2^{100}`` vertices. Computing the support function of
`M*X` for any square matrix `M` along the canonical direction
`e_{50} = [0, …, 0, 1, 0, …, 0]` takes around `10us` on a modern computer.

```@example tour
using LazySets.Arrays: SingleEntryVector

using BenchmarkTools

d = SingleEntryVector(50, 100, 1.0)
M = rand(100, 100)  # a random 100×100 matrix
H = rand(Hyperrectangle, dim=100)  # a random 100-dimensional hyperrectangle

out = @benchmark ρ($d, $M * $H)
out
```


## The Lazy paradigm

For the several set representations, such as polyhedra in constraint and in
vertex representation, ellipsoids, balls in different norms, and specific
classes of polyhedra (such as intervals, zonotopes, or hyperrectangles),
LazySets implements specialized methods for many set operations. Applications
that use LazySets can explore different approaches with minimal changes in their
code. Conversion between set representations as well as overapproximation and
underapproximation functionality are available.

Set operations can be performed in two possible (complementary) ways:
concretely or lazily. A *concrete* operation returns a set of a dedicated type
(such as `Hyperrectangle`). A *lazy* operation simply returns a wrapper object
representing the result of the operation between the given sets. For instance,
consider the transformation with a linear map, i.e., the set
``Y = \{y : y = Ax \textrm{ for some } x ∈ X\}``. In LazySets,
`linear_map(A, X)` returns a *concrete* set representation of ``Y``.
The algorithm that is actually used, as well as the type of ``Y``, depend on the
types of its arguments (this is the *multiple dispatch* paradigm) and the
dimension of `X`.

For example, if `X` is a polygon in vertex representation (`VPolygon`),
then `linear_map(M, X)` applies `M` to each vertex of `X`. However, if `X` is a
30-dimensional polyhedron in half-space representation (`HPolyhedron`), then the
half-space representation is used if `M` is invertible -- if not, and if `MX`
does not fall into any special case, the vertex representation has to be
computed. For polyhedra manipulation we use the library
[Polyhedra.jl](https://juliapolyhedra.github.io/Polyhedra.jl/dev/)
(the cost of converting between representations increases exponentially with
the dimension, so computing with concrete polyhedra is usually avoided in high
dimensions).
Observe that if we are interested in `P(MX)` for some projection matrix `P`,
that operation *can* be performed efficiently using the support function, i.e.,
without actually computing `MX`. On the other hand, there are special cases of
polyhedra for which concrete operations can be performed efficiently. For
instance, if `X` is a hyperrectangle (or a zonotope), the resulting set is a
zonotope and the computation is efficient.

To give an example of lazy operations, `LinearMap(A, X)`, or simply `A * X`,
computes the lazy linear map, which is just a new object that wraps the
computation of the linear map until it is actually needed. In other words,
`LinearMap(A, X)` can be used to reason about the linear map even if computing
the result is expensive (e.g., if `X` is high-dimensional), since this command
just builds an object representing the linear map of `A` and `X`. In LazySets,
Unicode symbols such as `A * X`, `X ⊕ Y`, `X ⊖ Y`, `X × Y`, all default to lazy
operations by design. Unicode is written using the LaTeX macro of the symbol
followed by the TAB key, such as `\oplus<tab>` for the (lazy) Minkowski sum ⊕.

Finally, one can combine lazy set operations to build lazy expressions that
represent several operations between sets, such as `Q = (Z ⊕ A*X) × T`.
By means of the basic tools of (convex) geometry, useful information about `Q`
can be obtained without actually computing the linear map, Minkowski sum and
Cartesian product in the above computation.
For example, if the computation only involves querying information about `Q` in
a restricted number of directions, as it is often the case in applications, the
lazy approach can be realized very efficiently!


## How to contribute

Development happens on [github](https://github.com/JuliaReach/LazySets.jl/).
New contributors should follow the links provided in the
[About](https://juliareach.github.io/LazySets.jl/dev/about/) section of this
documentation.


## How to cite

When citing `LazySets.jl`, please use the entry in [CITATION.bib](https://github.com/JuliaReach/LazySets.jl/blob/master/CITATION.bib).

Further publications using LazySets can be found in the [Publications](https://github.com/JuliaReach/LazySets.jl#-publications) section of the README. If you would like to list your work, feel free to create a [pull request](https://github.com/JuliaReach/LazySets.jl/pulls).
