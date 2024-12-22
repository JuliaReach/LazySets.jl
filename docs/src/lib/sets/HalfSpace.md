```@meta
CurrentModule = LazySets.HalfSpaceModule
```

# [Half-space (HalfSpace)](@id def_HalfSpace)

```@docs
HalfSpace
LinearConstraint
```

## Conversion

```julia
convert(::Type{HalfSpace{N,Vector{N}}}, hs::HalfSpace{N,<:AbstractVector{N}}) where {N}
```

The following method requires the [`SymEngine`](https://github.com/symengine/SymEngine.jl) package.

```@docs
LazySets.convert(::Type{HalfSpace{N}}, ::Expr; vars::Vector{Basic}=Basic[]) where {N}
```

## Operations

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
an_element(::LazySet)
```
```@meta
CurrentModule = LazySets.HalfSpaceModule
```
```@docs
an_element(::HalfSpace)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
complement(::LazySet)
```
```@meta
CurrentModule = LazySets.HalfSpaceModule
```
```@docs
complement(::HalfSpace)
constrained_dimensions(::HalfSpace)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isuniversal(::LazySet, ::Bool=false)
```
```@meta
CurrentModule = LazySets.HalfSpaceModule
```
```@docs
isuniversal(::HalfSpace, ::Bool=false)
normalize(::HalfSpace{N}, p::Real=N(2)) where {N}
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
rand(::Type{LazySet})
```
```@meta
CurrentModule = LazySets.HalfSpaceModule
```
```@docs
rand(::Type{HalfSpace})
distance(::AbstractVector, ::HalfSpace)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
∈(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.HalfSpaceModule
```
```@docs
∈(::AbstractVector, ::HalfSpace)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
project(::LazySet, ::AbstractVector{Int})
```
```@meta
CurrentModule = LazySets.HalfSpaceModule
```
```@docs
project(::HalfSpace, ::AbstractVector{Int})
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
ρ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.HalfSpaceModule
```
```@docs
ρ(::AbstractVector, ::HalfSpace)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
σ(::AbstractVector, ::LazySet)
```
```@meta
CurrentModule = LazySets.HalfSpaceModule
```
```@docs
σ(::AbstractVector, ::HalfSpace)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
translate(::LazySet, ::AbstractVector)
```
```@meta
CurrentModule = LazySets.HalfSpaceModule
```
```@docs
translate(::HalfSpace, ::AbstractVector)
iscomplement(::HalfSpace, ::HalfSpace)
```
```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
isdisjoint(::LazySet, ::LazySet)
```
```@meta
CurrentModule = LazySets.HalfSpaceModule
```
```@docs
isdisjoint(::HalfSpace, ::HalfSpace)
```

```@docs
halfspace_left(::AbstractVector, ::AbstractVector)
halfspace_right(::AbstractVector, ::AbstractVector)
HalfSpaceModule._ishalfspace
```

```@meta
CurrentModule = LazySets.API
```

Undocumented implementations:

* [`constraints_list`](@ref constraints_list(::LazySet))
* [`dim`](@ref dim(::LazySet))
* [`isbounded`](@ref isbounded(::LazySet))
* [`isempty`](@ref isempty(::LazySet))
* [`isoperationtype`](@ref isoperationtype(::Type{LazySet}))
* [`permute`](@ref permute(::LazySet, ::AbstractVector{Int}))

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`area`](@ref area(::LazySet))
* [`concretize`](@ref concretize(::LazySet))
* [`constraints`](@ref constraints(::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet))
* `copy(::Type{LazySet})`
* [`diameter`](@ref diameter(::LazySet, ::Real=Inf))
* [`eltype`](@ref eltype(::Type{<:LazySet}))
* [`eltype`](@ref eltype(::LazySet))
* [`isboundedtype`](@ref isboundedtype(::Type{LazySet}))
* [`isoperation`](@ref isoperation(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real=Inf))
* [`radius`](@ref radius(::LazySet, ::Real=Inf))
* [`rectify`](@ref rectify(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`surface`](@ref surface(::LazySet))
* [`vertices`](@ref vertices(::LazySet))
* [`affine_map`](@ref affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector))
* [`exponential_map`](@ref exponential_map(::AbstractMatrix, ::LazySet))
* [`is_interior_point`](@ref is_interior_point(::AbstractVector, ::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet))
* [`sample`](@ref sample(::LazySet, ::Int=1))
* [`scale`](@ref scale(::Real, ::LazySet))
* [`cartesian_product`](@ref cartesian_product(::LazySet, ::LazySet))
* [`convex_hull`](@ref convex_hull(::LazySet, ::LazySet))
* [`difference`](@ref difference(::LazySet, ::LazySet))
* [`distance`](@ref distance(::LazySet, ::LazySet; ::Real=2.0))
* [`exact_sum`](@ref exact_sum(::LazySet, ::LazySet))
* [`≈`](@ref ≈(::LazySet, ::LazySet))
* [`==`](@ref ==(::LazySet, ::LazySet))
* [`isequivalent`](@ref isequivalent(::LazySet, ::LazySet))
* [`⊂`](@ref ⊂(::LazySet, ::LazySet))
* [`minkowski_difference`](@ref minkowski_difference(::LazySet, ::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`linear_combination`](@ref linear_combination(::ConvexSet, ::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`extrema`](@ref extrema(::AbstractPolyhedron))
* [`extrema`](@ref extrema(::AbstractPolyhedron, ::Int))
* [`high`](@ref high(::AbstractPolyhedron))
* [`high`](@ref high(::AbstractPolyhedron, ::Int))
* [`isconvextype`](@ref isconvextype(::Type{AbstractPolyhedron}))
* [`ispolyhedral`](@ref ispolyhedral(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron))
* [`low`](@ref low(::AbstractPolyhedron, ::Int))
* [`vertices_list`](@ref vertices_list(::AbstractPolyhedron))
* [`intersection`](@ref intersection(::AbstractPolyhedron, ::AbstractPolyhedron))
* [`⊆`](@ref ⊆(::LazySet, ::AbstractPolyhedron))
* [`minkowski_sum`](@ref minkowski_sum(::AbstractPolyhedron, ::AbstractPolyhedron))
