# Set Interfaces

This section of the manual describes the interfaces for different set types.
Every set that fits the description of an interface should also implement it.
This helps in several ways:
- avoid code duplicates,
- provide functions for many sets at once,
- allow changes in the source code without changing the API.

The interface functions are outlined in the interface documentation.
For implementations of the interfaces see the corresponding sub-pages linked in
the respective sections.

!!! note
    The naming convention is such that all interface names (with the exception
    of the main abstract type `LazySet`) should be preceded by `Abstract`.

The following diagram shows the interface hierarchy.

![../assets/interfaces.png](../assets/interfaces.png)

```@contents
Pages = ["interfaces.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

## [General sets (LazySet)](@id def_LazySet)

Every convex set in this library implements this interface.

```@docs
LazySet
```

### Support function and support vector

Every `LazySet` type must define a function `σ` to compute the support vector.

```@docs
support_vector
ρ(::AbstractVector{N}, ::LazySet{N}) where {N<:Real}
support_function
σ
```

### Other globally defined set functions

```@docs
basetype
norm(::LazySet, ::Real=Inf)
radius(::LazySet, ::Real=Inf)
diameter(::LazySet, ::Real=Inf)
isbounded(::LazySet)
isbounded_unit_dimensions(::LazySet{N}) where {N<:Real}
an_element(::LazySet{N}) where {N<:Real}
tosimplehrep(::LazySet)
isuniversal(::LazySet{N}, ::Bool=false) where {N<:Real}
affine_map(M::AbstractMatrix, X::LazySet, v::AbstractVector)
reflect(::LazySet)
is_interior_point(::AbstractVector{N}, ::LazySet{N}; p=Inf, ε=_rtol(N)) where {N<:Real}
isoperationtype(::Type{<:LazySet})
isoperation(::LazySet)
isequivalent(::LazySet, ::LazySet)
isconvextype(::Type{<:LazySet})
```

Plotting is available for general one- or two-dimensional `LazySet`s, provided
that the overapproximation using iterative refinement is available:

```@docs
plot_recipe(::LazySet{N}, ::N=N(1e-3)) where {N<:Real}
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::LazySet{N}, ::N=N(1e-3)) where {N<:Real}
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::AbstractVector{VN}, ::N=N(1e-3), ::Int=40, ::Bool=false) where {N<:Real, VN<:LazySet{N}}
```

For three-dimensional sets, we support `Makie`:

```@docs
plot3d
plot3d!
```

### Set functions that override Base functions

```@docs
==(::LazySet, ::LazySet)
≈(::LazySet, ::LazySet)
copy(::LazySet)
eltype
```

### Aliases for set types

```@docs
CompactSet
NonCompactSet
```

### Implementations

Concrete set representations:

* [Empty set (EmptySet)](@ref def_EmptySet)
* [Polynomial zonotope (PolynomialZonotope)](@ref def_PolynomialZonotope)

Lazy set operations:

* [Affine map (AffineMap)](@ref def_AffineMap)
* [Linear map (LinearMap)](@ref def_LinearMap)
* [Exponential map (ExponentialMap)](@ref def_ExponentialMap)
* [Exponential projection map (ExponentialProjectionMap)](@ref def_ExponentialProjectionMap)
* [Reset map (ResetMap)](@ref def_ResetMap)
* [Translation](@ref def_Translation)
* [Bloating](@ref def_Bloating)
* [Binary Cartesian product (CartesianProduct)](@ref def_CartesianProduct)
* [``n``-ary Cartesian product (CartesianProductArray)](@ref def_CartesianProductArray)
* [Binary convex hull (ConvexHull)](@ref def_ConvexHull)
* [``n``-ary convex hull (ConvexHullArray)](@ref def_ConvexHullArray)
* [Binary intersection](@ref def_Intersection)
* [``n``-ary intersection (IntersectionArray)](@ref def_IntersectionArray)
* [Binary Minkowski sum (MinkowskiSum)](@ref def_MinkowskiSum)
* [``n``-ary Minkowski sum (MinkowskiSumArray)](@ref def_MinkowskiSumArray)
* [``n``-ary Minkowski sum with cache (CachedMinkowskiSumArray)](@ref def_CachedMinkowskiSumArray)
* [Binary set union (UnionSet)](@ref def_UnionSet)
* [``n``-ary set union (UnionSetArray)](@ref def_UnionSetArray)
* [Complement](@ref def_Complement)
* [Rectification](@ref def_Rectification)

## [Centrally symmetric sets (AbstractCentrallySymmetric)](@id def_AbstractCentrallySymmetric)

Centrally symmetric sets such as balls of different norms are characterized by a
center.
Note that there is a special interface combination
[Centrally symmetric polytope](@ref def_AbstractCentrallySymmetricPolytope).

```@docs
AbstractCentrallySymmetric
```

This interface defines the following functions:

```@docs
dim(::AbstractCentrallySymmetric)
isbounded(::AbstractCentrallySymmetric)
isuniversal(::AbstractCentrallySymmetric{N}, ::Bool=false) where {N<:Real}
an_element(::AbstractCentrallySymmetric{N}) where {N<:Real}
isempty(::AbstractCentrallySymmetric)
```

### Implementations

* [Euclidean-norm ball (Ball2)](@ref def_Ball2)
* [Ellipsoid](@ref def_Ellipsoid)
* [p-norm ball (Ballp)](@ref def_Ballp)

## [Polyhedra (AbstractPolyhedron)](@id def_AbstractPolyhedron)

A polyhedron has finitely many facets (*H-representation*) and is not
necessarily bounded.

```@docs
AbstractPolyhedron
```

This interface defines the following functions:

```@docs
∈(::AbstractVector{N}, ::AbstractPolyhedron{N}) where {N<:Real}
isuniversal(::AbstractPolyhedron{N}, ::Bool=false) where {N<:Real}
constrained_dimensions(::AbstractPolyhedron)
linear_map(::AbstractMatrix{N}, ::AbstractPolyhedron{N}) where {N<:Real}
chebyshev_center(::AbstractPolyhedron{N}) where {N<:AbstractFloat}
an_element(::AbstractPolyhedron{N}) where {N<:Real}
```

Plotting (bounded) polyhedra is available, too:

```@docs
plot_recipe(::AbstractPolyhedron{N}, ::N=zero(N)) where {N<:Real}
```

### Implementations

* [Half-space (HalfSpace)](@ref def_HalfSpace)
* [Polyhedron in constraint representation (HPolyhedron)](@ref def_HPolyhedron)
* [Hyperplane](@ref def_Hyperplane)
* [Line](@ref def_Line)
* [Universe](@ref def_Universe)

## [Polytopes (AbstractPolytope)](@id def_AbstractPolytope)

A polytope is a bounded set with finitely many vertices (*V-representation*)
resp. facets (*H-representation*).
Note that there is a special interface combination
[Centrally symmetric polytope](@ref def_AbstractCentrallySymmetricPolytope).

```@docs
AbstractPolytope
```

This interface defines the following functions:

```@docs
isbounded(::AbstractPolytope)
isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N<:Real}
singleton_list(::AbstractPolytope{N}) where {N<:Real}
isempty(::AbstractPolytope)
```

### Implementations

* [Polytope in constraint representation (HPolytope)](@ref def_HPolytope)
* [Polytope in vertex representation (VPolytope)](@ref def_VPolytope)

## [Polygons (AbstractPolygon)](@id def_AbstractPolygon)

A polygon is a two-dimensional polytope.

```@docs
AbstractPolygon
```

This interface defines the following functions:

```@docs
dim(P::AbstractPolygon)
linear_map(::AbstractMatrix{N}, P::AbstractPolygon{N}) where {N<:Real}
```

The following helper functions are used for sorting directions:

```@docs
LazySets.jump2pi
<=(::AbstractVector{N}, ::AbstractVector{N}) where {N<:AbstractFloat}
<=(::AbstractVector{N}, ::AbstractVector{N}) where {N<:Real}
LazySets.quadrant(::AbstractVector{Real})
```

### Implementations

* [Polygon in vertex representation (VPolygon)](@ref def_VPolygon)

## [Polygons in constraint representation (AbstractHPolygon)](@id def_AbstractHPolygon)

An HPolygon is a polygon in H-representation (or constraint representation).

```@docs
AbstractHPolygon
```

This interface defines the following functions:

```@docs
an_element(::AbstractHPolygon{N}) where {N<:Real}
∈(::AbstractVector{N}, ::AbstractHPolygon{N}) where {N<:Real}
rand(::Type{HPOLYGON}) where {HPOLYGON<:AbstractHPolygon}
tohrep(::HPOLYGON) where {HPOLYGON<:AbstractHPolygon}
tovrep(::AbstractHPolygon{N}) where {N<:Real}
addconstraint!(::AbstractHPolygon{N}, ::LinearConstraint{N}) where {N<:Real}
addconstraint!(::Vector{LC}, ::LinearConstraint{N}) where {N<:Real, LC<:LinearConstraint{N}}
isredundant(::LinearConstraint{N}, ::LinearConstraint{N}, ::LinearConstraint{N}) where {N<:Real}
remove_redundant_constraints!(::AbstractHPolygon)
constraints_list(::AbstractHPolygon{N}) where {N<:Real}
vertices_list(::AbstractHPolygon{N}) where {N<:Real}
isbounded(::AbstractHPolygon, ::Bool=true)
```

### Implementations

* [Polygon in constraint representation (HPolygon)](@ref def_HPolygon)
* [Polygon in optimized constraint representation (HPolygonOpt)](@ref def_HPolygonOpt)

## [Centrally symmetric polytopes (AbstractCentrallySymmetricPolytope)](@id def_AbstractCentrallySymmetricPolytope)

A centrally symmetric polytope is a combination of two other interfaces:
[Centrally symmetric sets](@ref def_AbstractCentrallySymmetric) and
[Polytope](@ref def_AbstractPolytope).

```@docs
AbstractCentrallySymmetricPolytope
```

This interface defines the following functions:

```@docs
dim(::AbstractCentrallySymmetricPolytope)
an_element(::AbstractCentrallySymmetricPolytope{N}) where {N<:Real}
isempty(::AbstractCentrallySymmetricPolytope)
```

### Implementations

* [Manhattan-norm ball (Ball1)](@ref def_Ball1)

## [Zonotopes (AbstractZonotope)](@id def_AbstractZonotope)

A zonotope is a specific centrally symmetric polytope characterized by a
center and a collection of generators.

```@docs
AbstractZonotope
```

This interface defines the following functions:

```@docs
ngens(::AbstractZonotope)
genmat_fallback(::AbstractZonotope{N}) where {N<:Real}
generators_fallback(::AbstractZonotope{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::AbstractZonotope{N}) where {N<:Real}
σ(::AbstractVector{N}, ::AbstractZonotope{N}) where {N<:Real}
∈(::AbstractVector{N}, ::AbstractZonotope{N}) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::AbstractZonotope{N}) where {N<:Real}
translate(::AbstractZonotope{N}, ::AbstractVector{N}) where {N<:Real}
constraints_list(::AbstractZonotope{N}) where {N<:Real}
constraints_list(::AbstractZonotope{N}; ::Bool=true) where {N<:AbstractFloat}
vertices_list(::AbstractZonotope{N}) where {N<:Real}
order(::AbstractZonotope)
```

### Implementations

* [Zonotope](@ref def_Zonotope)
* [Line segment (LineSegment)](@ref def_LineSegment)

## [Hyperrectangles (AbstractHyperrectangle)](@id def_AbstractHyperrectangle)

A hyperrectangle is a special centrally symmetric polytope with axis-aligned
facets.

```@docs
AbstractHyperrectangle
```

This interface defines the following functions:

```@docs
norm(::AbstractHyperrectangle, ::Real=Inf)
radius(::AbstractHyperrectangle, ::Real=Inf)
σ(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real}
∈(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real}
vertices_list(::AbstractHyperrectangle{N}) where {N<:Real}
constraints_list(::AbstractHyperrectangle{N}) where {N<:Real}
high(::AbstractHyperrectangle{N}) where {N<:Real}
high(::AbstractHyperrectangle{N}, ::Int) where {N<:Real}
low(::AbstractHyperrectangle{N}) where {N<:Real}
low(::AbstractHyperrectangle{N}, ::Int) where {N<:Real}
isflat(::AbstractHyperrectangle)
split(::AbstractHyperrectangle{N}, ::AbstractVector{Int}) where {N<:Real}
generators(::AbstractHyperrectangle)
genmat(::AbstractHyperrectangle)
ngens(::AbstractHyperrectangle{N}) where {N<:Real}
rectify(::AbstractHyperrectangle)
```

### Implementations

Concrete set representations:

* [Hyperrectangle](@ref def_Hyperrectangle)
* [Infinity-norm ball (BallInf)](@ref def_BallInf)
* [Interval](@ref def_Interval)

Lazy set operations:

* [Symmetric interval hull (SymmetricIntervalHull)](@ref def_SymmetricIntervalHull)

## [Singletons (AbstractSingleton)](@id def_AbstractSingleton)

A singleton is a special hyperrectangle consisting of only one point.

```@docs
AbstractSingleton
```

This interface defines the following functions:

```@docs
σ(::AbstractVector{N}, ::AbstractSingleton{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::AbstractSingleton{N}) where {N<:Real}
∈(::AbstractVector{N}, ::AbstractSingleton{N}) where {N<:Real}
an_element(::AbstractSingleton{N}) where {N<:Real}
center(::AbstractSingleton{N}) where {N<:Real}
vertices_list(::AbstractSingleton{N}) where {N<:Real}
radius_hyperrectangle(::AbstractSingleton{N}) where {N<:Real}
radius_hyperrectangle(::AbstractSingleton{N}, ::Int) where {N<:Real}
high(::AbstractSingleton{N}) where {N<:Real}
high(::AbstractSingleton{N}, ::Int) where {N<:Real}
low(::AbstractSingleton{N}) where {N<:Real}
low(::AbstractSingleton{N}, ::Int) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::AbstractSingleton{N}) where {N<:Real}
generators(::AbstractSingleton{N}) where {N<:Real}
genmat(::AbstractSingleton{N}) where {N<:Real}
ngens(::AbstractSingleton)
plot_recipe(::AbstractSingleton{N}, ::N=zero(N)) where {N<:Real}
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::AbstractSingleton{N}, ::N=zero(N)) where {N<:Real}
```

### Implementations

* [Singleton](@ref def_Singleton)
* [Origin (ZeroSet)](@ref def_ZeroSet)
