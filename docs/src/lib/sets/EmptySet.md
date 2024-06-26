```@meta
CurrentModule = LazySets.EmptySetModule
```

# [Empty set (EmptySet)](@id def_EmptySet)

```@docs
EmptySet
∅
```

## Operations

```@docs
an_element(::EmptySet)
area(::EmptySet)
chebyshev_center_radius(::EmptySet; kwargs...)
complement(::EmptySet{N}) where {N}
diameter(::EmptySet, ::Real=Inf)
dim(::EmptySet)
high(::EmptySet)
high(::EmptySet, ::Int)
isbounded(::EmptySet)
isempty(::EmptySet)
isuniversal(::EmptySet{N}, ::Bool=false) where {N}
low(::EmptySet)
low(::EmptySet, ::Int)
norm(::EmptySet, ::Real=Inf)
radius(::EmptySet, ::Real=Inf)
rand(::Type{EmptySet})
rectify(::EmptySet)
reflect(::EmptySet)
vertices_list(::EmptySet)
vertices(::EmptySet)
volume(::EmptySet{N}) where {N}
∈(::AbstractVector, ::EmptySet)
linear_map(::AbstractMatrix, ::EmptySet)
ρ(::AbstractVector, ::EmptySet)
σ(::AbstractVector, ::EmptySet)
translate(::EmptySet, ::AbstractVector)
plot_recipe(::EmptySet{N}, ::Any=zero(N)) where {N}
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:
* [`convex_hull`](@ref convex_hull(::LazySet))
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`isconvextype`](@ref isconvextype(::Type{LazySet}))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`scale!`](@ref scale!(::Real, ::LazySet))
* [`translate!`](@ref translate!(::LazySet, ::AbstractVector))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`intersection`](@ref isequivalent(::LazySet, ::LazySet))
* [`isdisjoint`](@ref isdisjoint(::LazySet, ::LazySet))
* [`⊆`](@ref ⊆(::LazySet, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`concretize`](@ref concretize(::LazySet))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`extrema`](@ref extrema(::LazySet))
* [`extrema`](@ref extrema(::LazySet, ::Int))
* [`isoperation`](@ref isoperation(::LazySet))
* [`is_polyhedral`](@ref is_polyhedral(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
