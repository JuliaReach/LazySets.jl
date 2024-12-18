```@meta
CurrentModule = LazySets.UniverseModule
```

# [Universe](@id def_Universe)

```@docs
Universe
```

## Operations

```@docs
an_element(::Universe{N}) where {N}
complement(::Universe{N}) where {N}
constrained_dimensions(::Universe)
constraints(::Universe{N}) where {N}
constraints_list(::Universe{N}) where {N}
diameter(::Universe, ::Real=Inf)
dim(::Universe)
isbounded(::Universe)
isempty(::Universe)
isuniversal(::Universe{N}, ::Bool=false) where {N}
norm(::Universe, ::Real=Inf)
radius(::Universe, ::Real=Inf)
rand(::Type{Universe})
reflect(::Universe)
∈(::AbstractVector, ::Universe)
permute(::Universe, ::AbstractVector{Int})
ρ(::AbstractVector, ::Universe)
σ(::AbstractVector, ::Universe)
translate(::Universe, ::AbstractVector)
translate!(::Universe, ::AbstractVector)
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* `copy(::Universe)`
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`project`](@ref project(::LazySet, ::AbstractVector))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`scale!`](@ref scale!(::Real, ::LazySet))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`intersection`](@ref intersection(::LazySet, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`AbstractPolyhedron`](@ref):
* [`extrema`](@ref extrema(::LazySet))
* [`extrema`](@ref extrema(::LazySet, ::Int))
* [`high`](@ref high(::AbstractPolyhedron))
* [`high`](@ref high(::AbstractPolyhedron, ::Int))
* [`isconvextype`](@ref isconvextype(::Type{<:AbstractPolyhedron}))
* [`ispolyhedral`](@ref ispolyhedral(::LazySet))
* [`low`](@ref low(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron, ::Int))
* [`vertices_list`](@ref vertices_list(::AbstractPolyhedron))
* [`vertices`](@ref vertices(::AbstractPolyhedron))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractPolyhedron, ::AbstractPolyhedron))

Inherited from [`ConvexSet`](@ref):
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`LazySet`](@ref):
* [`area`](@ref area(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isoperation`](@ref isoperation(::LazySet))
* [`rectify`](@ref rectify(::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`isdisjoint`](@ref isdisjoint(::LazySet, ::LazySet))
* [`⊆`](@ref ⊆(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
