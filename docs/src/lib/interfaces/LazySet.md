```@contents
Pages = ["LazySet.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [General sets (LazySet)](@id def_LazySet)

Every set type in this library is a subtype of the abstract type `LazySet`.

```@docs
LazySet
```

## Plotting

Plotting via the `Plots` package is available for one- or two-dimensional sets.
The default algorithm is to plot an outer approximation using the support
function (1D) respectively the support vector (2D). This means that (1) plotting
will fail if these functionalities are not available (e.g., for lazy
`Intersection`s) and (2) that plots of non-convex sets can be misleading. The
implementation below internally relies on the function `plot_recipe`. For some
set types (e.g., `Intersection`), the default implementation is overridden.

```@docs
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::LazySet{N}, ::Real=N(1e-3)) where {N}
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::AbstractVector{VN}, ::Real=N(1e-3), ::Int=40; ::Bool=false) where {N, VN<:LazySet{N}}
plot_vlist(::S, ::Real) where {S<:LazySet}
```

For three-dimensional sets, we support `Makie`:

```@docs
plot3d
plot3d!
```

## Globally defined set functions

```@docs
○(c, a)
isconvextype(::Type{<:LazySet})
low(::LazySet)
high(::LazySet)
extrema(::LazySet, ::Int)
extrema(::LazySet)
convex_hull(::LazySet; kwargs...)
triangulate(::LazySet)
isboundedtype(::Type{<:LazySet})
isbounded(::LazySet)
_isbounded_unit_dimensions(::LazySet)
is_polyhedral(::LazySet)
isfeasible
norm(::LazySet, ::Real=Inf)
radius(::LazySet, ::Real=Inf)
diameter(::LazySet, ::Real=Inf)
isempty(::LazySet{N}, ::Bool=false) where {N}
linear_map(::AbstractMatrix, ::LazySet; kwargs...)
linear_map(::Number, ::LazySet; kwargs...)
affine_map(::Any, ::LazySet, ::AbstractVector)
exponential_map(::AbstractMatrix, ::LazySet)
an_element(::LazySet)
tosimplehrep(::LazySet)
reflect(::LazySet)
is_interior_point(::AbstractVector{<:Real}, ::LazySet; kwargs...)
isoperation(::LazySet)
isequivalent(::LazySet, ::LazySet)
surface(::LazySet)
area(::LazySet)
concretize(::LazySet)
complement(::LazySet)
polyhedron(::LazySet)
project(::LazySet, ::AbstractVector{Int}, ::Nothing=nothing, ::Int=dim(S))
project(::LazySet, ::AbstractVector{Int}, ::Type{TS}, ::Int=dim(S)) where {TS<:LazySet}
project(::LazySet, ::AbstractVector{Int}, ::Pair{T, N}, ::Int=dim(S)) where {T<:UnionAll, N<:Real}
project(::LazySet, ::AbstractVector{Int}, ::Real, ::Int=dim(S))
rectify(::LazySet, ::Bool=false)
permute
rationalize(::Type{T}, ::LazySet{<:AbstractFloat}, ::Real) where {T<:Integer}
singleton_list(::LazySet)
constraints(::LazySet)
vertices(::LazySet)
delaunay
chebyshev_center_radius(::LazySet{N}) where {N}
scale(::Real, ::LazySet)
plot_recipe(::LazySet, ::Any)
```

The following methods are also defined for `LazySet` but cannot be documented
due to a bug in the documentation package.

```@docs
low(::ConvexSet{N}, ::Int) where {N}
high(::ConvexSet{N}, ::Int) where {N}
an_element(::ConvexSet{N}) where {N}
```

## Support function and support vector

Every `LazySet` type must define a function `σ` to compute the support vector.
The support function, `ρ`, can optionally be defined; otherwise, a fallback
definition based on `σ` is used.

```@docs
σ
ρ(::AbstractVector, ::LazySet)
```

## Set functions that override Base functions

```@docs
==(::LazySet, ::LazySet)
≈(::LazySet, ::LazySet)
```

## Aliases for set types

```@docs
CompactSet
NonCompactSet
```

## Implementations

Concrete set representations:

* [Empty set (EmptySet)](@ref def_EmptySet)

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
