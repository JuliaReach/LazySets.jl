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

This interface requires to implement the following function:

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
dim(::LazySet)
```
```@meta
CurrentModule = LazySets
```

## Support vector and support function

Most `LazySet` types (particularly convex ones) define a function `σ` to compute
the support vector.
The support function, `ρ`, can optionally be defined; otherwise, a fallback
definition based on `σ` is used.

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
σ(::AbstractVector, ::LazySet)
ρ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
ρ(::AbstractVector, ::LazySet)
```

## Globally defined set functions and default implementations

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
an_element(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
an_element(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
area(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
area(::LazySet)
chebyshev_center_radius(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
complement(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
complement(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
concretize(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
concretize(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
constraints(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
constraints(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
convex_hull(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
convex_hull(::LazySet; kwargs...)
```
`copy(::Type{LazySet})`
```@docs
triangulate(X::LazySet; ::String="delaunay")
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
diameter(::LazySet, ::Real=Inf)
```
```@meta
CurrentModule = LazySets
```
```@docs
diameter(::LazySet, ::Real=Inf)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
eltype(::Type{<:LazySet})
```
```@meta
CurrentModule = LazySets
```
```@docs
eltype(::Type{<:LazySet})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
eltype(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
eltype(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
extrema(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
extrema(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
extrema(::LazySet, ::Int)
```
```@meta
CurrentModule = LazySets
```
```@docs
extrema(::LazySet, ::Int)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
high(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
high(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isbounded(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
isbounded(::LazySet)
_isbounded_unit_dimensions(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isboundedtype(::Type{<:LazySet})
```
```@meta
CurrentModule = LazySets
```
```@docs
isboundedtype(::Type{<:LazySet})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isconvextype(::Type{<:LazySet})
```
```@meta
CurrentModule = LazySets
```
```@docs
isconvextype(::Type{<:LazySet})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isempty(::LazySet, ::Bool=false)
```
```@meta
CurrentModule = LazySets
```
```@docs
isempty(::LazySet, ::Bool=false)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isoperation(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
isoperation(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
ispolyhedral(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
ispolyhedral(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
ispolyhedraltype(::Type{<:LazySet})
```
```@meta
CurrentModule = LazySets
```
```@docs
ispolyhedraltype(::Type{<:LazySet})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
ispolytopic(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
ispolytopic(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
ispolytopictype(::Type{<:LazySet})
```
```@meta
CurrentModule = LazySets
```
```@docs
ispolytopictype(::Type{<:LazySet})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
low(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
low(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
norm(::LazySet, ::Real=Inf)
```
```@meta
CurrentModule = LazySets
```
```@docs
norm(::LazySet, ::Real=Inf)
polyhedron(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
radius(::LazySet, ::Real=Inf)
```
```@meta
CurrentModule = LazySets
```
```@docs
radius(::LazySet, ::Real=Inf)
rationalize(::Type{T}, ::LazySet{<:AbstractFloat}, ::Real) where {T<:Integer}
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rectify(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
rectify(::LazySet, ::Bool=false)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
reflect(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
reflect(::LazySet)
singleton_list(::LazySet)
triangulate_faces(::LazySet)
tohrep(::LazySet)
tosimplehrep(::LazySet)
tovrep(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
vertices(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
vertices(::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets
```
```@docs
affine_map(::Any, ::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
exponential_map(::AbstractMatrix, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
exponential_map(::AbstractMatrix, ::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
is_interior_point(::AbstractVector{<:Real}, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
is_interior_point(::AbstractVector{N}, ::LazySet{N}; p=N(Inf), ε=_rtol(N)) where {N<:Real}
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
linear_map(::AbstractMatrix, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
linear_map(::AbstractMatrix, ::LazySet; kwargs...)
project(::LazySet, ::AbstractVector{Int}, ::Nothing=nothing, ::Int=dim(X))
project(::LazySet, ::AbstractVector{Int}, ::Type{<:LazySet}, ::Int=dim(X))
project(::LazySet, ::AbstractVector{Int}, ::Pair{<:UnionAll,<:Real}, ::Int=dim(X))
project(::LazySet, ::AbstractVector{Int}, ::Real, ::Int=dim(X))
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
scale(::Real, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
scale(::Real, ::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
translate(::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets
```
```@docs
translate(::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
cartesian_product(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
cartesian_product(::LazySet, ::LazySet)
convex_hull(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isapprox(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
isapprox(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isdisjoint(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
isdisjoint(::LazySet, ::LazySet, ::Bool=false)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
==(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
==(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isequivalent(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
isequivalent(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
⊂(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
⊂(::LazySet, ::LazySet, ::Bool=false)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
issubset(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
issubset(::LazySet, ::LazySet, ::Bool=false)
minkowski_difference(::LazySet, ::LazySet)
minkowski_sum(::LazySet, ::LazySet)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`high`](@ref high(::LazySet, ::Int))
* [`low`](@ref low(::LazySet, ::Int))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`linear_combination`](@ref linear_combination(::LazySet, ::LazySet))

```@meta
CurrentModule = LazySets
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
LazySets.apply_recipe(::AbstractDict{Symbol,Any}, ::LazySet{N}, ::Real=N(1e-3)) where {N}
LazySets.apply_recipe(::AbstractDict{Symbol,Any}, ::AbstractVector{VN}, ::Real=N(1e-3), ::Int=40; ::Bool=false) where {N, VN<:LazySet{N}}
plot_vlist(::LazySet, ::Real)
plot_recipe(::LazySet, ::Any)
```

For three-dimensional sets, we support `Makie`:

```@docs
plot3d
plot3d!
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
