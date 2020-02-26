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
plot_recipe(::LineSegment{N}, ::N=zero(N)) where {N<:Real}
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::Union{LineSegment{N}, Interval{N}}, ::N=zero(N)) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))
