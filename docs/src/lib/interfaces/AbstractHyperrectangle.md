```@contents
Pages = ["AbstractHyperrectangle.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Hyperrectangles (AbstractHyperrectangle)](@id def_AbstractHyperrectangle)

A hyperrectangle is a special centrally symmetric polytope with axis-aligned
facets.

```@docs
AbstractHyperrectangle
```

This interface requires to implement the following function:

```@docs
radius_hyperrectangle(::AbstractHyperrectangle)
```

This interface defines the following functions:

```@docs
□(c, r)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
constraints_list(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
constraints_list(::AbstractHyperrectangle)
isflat(::AbstractHyperrectangle)
```
```@docs; canonical=false
ngens(::AbstractZonotope)
```
```@docs
ngens(::AbstractHyperrectangle)
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
norm(::AbstractHyperrectangle, ::Real=Inf)
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
radius(::AbstractHyperrectangle, ::Real=Inf)
radius_hyperrectangle(::AbstractHyperrectangle, ::Int)
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
rectify(::AbstractHyperrectangle)
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
reflect(::AbstractHyperrectangle)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
vertices_list(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
vertices_list(::AbstractHyperrectangle)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
volume(::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
volume(::AbstractHyperrectangle)
distance(::AbstractVector, ::AbstractHyperrectangle{N}; ::Real=N(2)) where {N}
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
∈(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
∈(::AbstractVector, ::AbstractHyperrectangle)
split(::AbstractHyperrectangle, ::AbstractVector{Int})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
σ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
σ(::AbstractVector, ::AbstractHyperrectangle)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
difference(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
difference(::AbstractHyperrectangle, ::AbstractHyperrectangle)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
intersection(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
intersection(::AbstractHyperrectangle, ::AbstractHyperrectangle)
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
isdisjoint(::AbstractHyperrectangle, ::AbstractHyperrectangle, ::Bool=false)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
⊆(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
⊆(::AbstractHyperrectangle, ::AbstractHyperrectangle, ::Bool=false)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
minkowski_difference(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
minkowski_difference(::AbstractHyperrectangle, ::AbstractHyperrectangle)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
minkowski_sum(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets
```
```@docs
minkowski_sum(::AbstractHyperrectangle, ::AbstractHyperrectangle)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`area`](@ref area(::LazySet))
* [`extrema`](@ref extrema(::LazySet))
* [`extrema`](@ref extrema(::LazySet, ::Int))
```@meta
CurrentModule = LazySets
```
* [`generators`](@ref generators(::AbstractZonotope))
* [`genmat`](@ref genmat(::AbstractZonotope))
```@meta
CurrentModule = LazySets.API
```
* [`high`](@ref high(::LazySet))
* [`high`](@ref high(::LazySet, ::Int))
* [`low`](@ref low(::LazySet))
* [`low`](@ref low(::LazySet, ::Int))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`distance`](@ref distance(::LazySet, ::LazySet; ::Real=2.0))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`chebyshev_center_radius`](@ref chebyshev_center_radius(::LazySet))
* [`complement`](@ref complement(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* [`constraints`](@ref constraints(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* `copy(::Type{LazySet})`
* [`delaunay`](@ref delaunay(::LazySet))
* [`diameter`](@ref diameter(::LazySet, ::Real=Inf))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`rationalize`](@ref rationalize(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`tosimplehrep`](@ref tosimplehrep(::LazySet))
* [`triangulate`](@ref triangulate(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{LazySet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isboundedtype`](@ref isboundedtype(::Type{AbstractPolytope}))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`center`](@ref center(::AbstractCentrallySymmetricPolytope, ::Int))
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope, ::Bool=false))

Inherited from [`AbstractZonotope`](@ref):
* [`order`](@ref order(::AbstractZonotope))
* [`remove_redundant_generators`](@ref remove_redundant_generators(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))
* [`reduce_order`](@ref reduce_order(::AbstractZonotope, ::Real, ::AbstractReductionMethod=GIR05()))
* [`split`](@ref split(::AbstractZonotope, ::Int))
* [`split`](@ref split(::AbstractZonotope, ::AbstractVector{Int}, ::AbstractVector{Int}))

## Implementations

Concrete set representations:

* [Hyperrectangle](@ref def_Hyperrectangle)
* [Infinity-norm ball (BallInf)](@ref def_BallInf)
* [Interval](@ref def_Interval)

Lazy set operations:

* [Symmetric interval hull (SymmetricIntervalHull)](@ref def_SymmetricIntervalHull)
