```@meta
CurrentModule = LazySets
```

# Intersection

## [Binary intersection (Intersection)](@id def_Intersection)

```@docs
Intersection
dim(::Intersection)
ρ(::AbstractVector, ::Intersection)
ρ(::AbstractVector, ::Intersection{N, S1, S2}) where {N, S1<:LazySet, S2<:Union{HalfSpace, Hyperplane, Line2D}}
ρ(::AbstractVector, ::Intersection{N, S1, S2}) where {N, S1<:LazySet, S2<:AbstractPolyhedron}
ρ(::AbstractVector, ::Intersection{N, S1, S2}) where {N, S1<:AbstractPolyhedron, S2<:AbstractPolyhedron}
σ(::AbstractVector, ::Intersection)
isbounded(::Intersection)
isempty(::Intersection)
in(::AbstractVector, ::Intersection)
constraints_list(::Intersection)
vertices_list(::Intersection)
isempty_known(::Intersection)
set_isempty!(::Intersection, ::Bool)
swap(::Intersection)
use_precise_ρ
_line_search
_projection
linear_map(::AbstractMatrix, ::Intersection)
plot_recipe(::Intersection{N}, ::N=zero(N), ::Int=40) where {N}
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::Intersection{N}, ::Real=zero(N), ::Int=40) where {N}
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet)
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`reflect`](@ref reflect(::LazySet))

### Intersection cache

```@docs
IntersectionCache
```

## [``n``-ary intersection (IntersectionArray)](@id def_IntersectionArray)

```@docs
IntersectionArray
Intersection!
dim(::IntersectionArray)
σ(::AbstractVector, ::IntersectionArray)
isbounded(::IntersectionArray)
in(::AbstractVector, ::IntersectionArray)
array(::IntersectionArray)
constraints_list(::IntersectionArray)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet)
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
