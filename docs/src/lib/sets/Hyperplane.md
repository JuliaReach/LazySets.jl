```@meta
CurrentModule = LazySets.HyperplaneModule
```

# [Hyperplane](@id def_Hyperplane)

```@docs
Hyperplane
```

## Conversion

The following method requires the [`SymEngine`](https://github.com/symengine/SymEngine.jl) package.

```@docs
LazySets.convert(::Type{Hyperplane{N}}, ::Expr; vars::Vector{Basic}=Basic[]) where {N}
```

## Operations

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
an_element(::LazySet)
```
```@meta
CurrentModule = LazySets.HyperplaneModule
```
```@docs
an_element(::Hyperplane)
constrained_dimensions(::Hyperplane)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isbounded(::LazySet)
```
```@meta
CurrentModule = LazySets.HyperplaneModule
```
```@docs
isbounded(::Hyperplane)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isuniversal(::LazySet)
```
```@meta
CurrentModule = LazySets.HyperplaneModule
```
```@docs
isuniversal(::Hyperplane, ::Bool=false)
normalize(::Hyperplane{N}, p::Real=N(2)) where {N}
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.HyperplaneModule
```
```@docs
rand(::Type{Hyperplane})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
∈(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.HyperplaneModule
```
```@docs
∈(::AbstractVector, ::Hyperplane)
project(::AbstractVector, ::Hyperplane)
reflect(::AbstractVector, ::Hyperplane)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
ρ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.HyperplaneModule
```
```@docs
ρ(::AbstractVector, ::Hyperplane)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
σ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.HyperplaneModule
```
```@docs
σ(::AbstractVector, ::Hyperplane)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
translate(::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.HyperplaneModule
```
```@docs
translate(::Hyperplane, ::AbstractVector)
```

```@docs
HyperplaneModule._ishyperplanar
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`constraints_list`](@ref constraints_list(::LazySet))
* [`dim`](@ref dim(::LazySet))
* [`isempty`](@ref isempty(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`distance`](@ref distance(::AbstractVector, ::LazySet))
* [`project`](@ref project(::LazySet, ::AbstractVector{Int}))
* [`isdisjoint`](@ref isdisjoint(::LazySet, ::LazySet))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`area`](@ref area(::LazySet))
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
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`isoperation`](@ref isoperation(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real=Inf))
* [`polyhedron`](@ref polyhedron(::LazySet))
* [`radius`](@ref radius(::LazySet, ::Real=Inf))
* [`rationalize`](@ref rationalize(::LazySet))
* [`rectify`](@ref rectify(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`tosimplehrep`](@ref tosimplehrep(::LazySet))
* [`triangulate`](@ref triangulate(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`isconvextype`](@ref isconvextype(::Type{ConvexSet}))
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`extrema`](@ref extrema(::AbstractPolyhedron))
* [`extrema`](@ref extrema(::AbstractPolyhedron, ::Int))
* [`high`](@ref high(::AbstractPolyhedron))
* [`high`](@ref high(::AbstractPolyhedron, ::Int))
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron, ::Int))
* [`vertices_list`](@ref vertices_list(::AbstractPolyhedron))
* [`intersection`](@ref intersection(::AbstractPolyhedron, ::AbstractPolyhedron))
* [`⊆`](@ref ⊆(::LazySet, ::AbstractPolyhedron))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractPolyhedron, ::AbstractPolyhedron))
