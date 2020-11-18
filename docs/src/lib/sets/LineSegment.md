```@meta
CurrentModule = LazySets
```

# [Line segment (LineSegment)](@id def_LineSegment)

```@docs
LineSegment
dim(::LineSegment)
σ(::AbstractVector{N}, ::LineSegment{N}) where {N<:Real}
∈(::AbstractVector{N}, ::LineSegment{N}) where {N<:Real}
center(::LineSegment{N}) where {N<:Real}
an_element(::LineSegment{N}) where {N<:Real}
rand(::Type{LineSegment})
halfspace_left(::LineSegment)
halfspace_right(::LineSegment)
vertices_list(::LineSegment{N}) where {N<:Real}
constraints_list(::LineSegment{N}) where {N<:Real}
translate(::LineSegment{N}, ::AbstractVector{N}) where {N<:Real}
generators(::LineSegment{N}) where {N<:Real}
genmat(::LineSegment)
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::Union{LineSegment{N}, Interval{N}}, ::N=zero(N)) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope, ::Bool=false))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))
