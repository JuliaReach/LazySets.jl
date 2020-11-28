```@meta
CurrentModule = LazySets
```

# Intersection

## [Binary intersection (Intersection)](@id def_Intersection)

```@docs
Intersection
∩(::LazySet, ::LazySet)
dim(::Intersection)
ρ(::AbstractVector, ::Intersection)
ρ(::AbstractVector, ::Intersection{N, S1, S2}) where {N, S1<:LazySet{N}, S2<:Union{HalfSpace{N}, Hyperplane{N}, Line2D{N}}}
ρ(::AbstractVector, ::Intersection{N, S1, S2}) where {N, S1<:LazySet{N}, S2<:AbstractPolyhedron{N}}
ρ(::AbstractVector, ::Intersection{N, S1, S2}) where {N, S1<:AbstractPolyhedron{N}, S2<:AbstractPolyhedron{N}}
σ(::AbstractVector, ::Intersection)
isbounded(::Intersection)
isempty(::Intersection)
∈(::AbstractVector, ::Intersection)
constraints_list(::Intersection)
isempty_known(::Intersection)
set_isempty!(::Intersection, ::Bool)
swap(::Intersection)
use_precise_ρ
_line_search
_projection
linear_map(::AbstractMatrix, ::Intersection)
plot_recipe(::Intersection{N}, ::N=zero(N), ::Int=40) where {N<:Real}
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::Intersection{N}, ::N=zero(N), ::Int=40) where {N}
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N})
* [`singleton_list`](@ref singleton_list(::LazySet))

### Intersection cache

```@docs
IntersectionCache
```

## [``n``-ary intersection (IntersectionArray)](@id def_IntersectionArray)

```@docs
IntersectionArray
dim(::IntersectionArray)
σ(::AbstractVector, ::IntersectionArray)
isbounded(::IntersectionArray)
∈(::AbstractVector, ::IntersectionArray)
array(::IntersectionArray)
constraints_list(::IntersectionArray{N}) where {N}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N})
* [`singleton_list`](@ref singleton_list(::LazySet))
