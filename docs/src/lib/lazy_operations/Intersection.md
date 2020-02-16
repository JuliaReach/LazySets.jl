```@meta
CurrentModule = LazySets
```

# Intersection

## [Binary intersection (Intersection)](@id def_Intersection)

```@docs
Intersection
∩(::LazySet, ::LazySet)
dim(::Intersection)
ρ(::AbstractVector{N}, ::Intersection{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::Intersection{N, S1, S2}) where {N<:Real, S1<:LazySet{N}, S2<:Union{HalfSpace{N}, Hyperplane{N}, Line{N}}}
ρ(::AbstractVector{N}, ::Intersection{N, S1, S2}) where {N<:Real, S1<:LazySet{N}, S2<:AbstractPolyhedron{N}}
ρ(::AbstractVector{N}, ::Intersection{N, S1, S2}) where {N<:Real, S1<:AbstractPolyhedron{N}, S2<:AbstractPolyhedron{N}}
σ(::AbstractVector{N}, ::Intersection{N}) where {N<:Real}
isbounded(::Intersection)
isempty(::Intersection)
∈(::AbstractVector{N}, ::Intersection{N}) where {N<:Real}
constraints_list(::Intersection{N}) where {N<:Real}
isempty_known(::Intersection)
set_isempty!(::Intersection, ::Bool)
swap(::Intersection{N, S1, S2}) where {N<:Real, S1, S2}
use_precise_ρ
_line_search
_projection
linear_map(::AbstractMatrix{N}, ::Intersection{N}) where {N}
plot_recipe(::Intersection{N}, ::N=zero(N), ::Int=40) where {N<:Real}
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::Intersection{N}, ::N=zero(N), ::Int=40) where {N<:Real}
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})

### Intersection cache

```@docs
IntersectionCache
```

## [``n``-ary intersection (IntersectionArray)](@id def_IntersectionArray)

```@docs
IntersectionArray
dim(::IntersectionArray)
σ(::AbstractVector{N}, ::IntersectionArray{N}) where {N<:Real}
isbounded(::IntersectionArray)
∈(::AbstractVector{N}, ::IntersectionArray{N}) where {N<:Real}
array(::IntersectionArray{N, S}) where {N<:Real, S<:LazySet{N}}
constraints_list(::IntersectionArray{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`an_element`](@ref an_element(::LazySet{N}) where {N<:Real})
