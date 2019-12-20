```@meta
CurrentModule = LazySets
```

# [Polygon in vertex representation (VPolygon)](@id def_VPolygon)

```@docs
VPolygon
σ(::AbstractVector{N}, ::VPolygon{N}) where {N<:Real}
∈(::AbstractVector{N}, ::VPolygon{N}) where {N<:Real}
an_element(::VPolygon{N}) where {N<:Real}
rand(::Type{VPolygon})
vertices_list(::VPolygon{N}) where {N<:Real}
tohrep(::VPolygon{N}, ::Type{HPOLYGON}=HPolygon) where {N<:Real, HPOLYGON<:AbstractHPolygon}
tovrep(::VPolygon{N}) where {N<:Real}
constraints_list(::VPolygon{N}) where {N<:Real}
translate(::VPolygon{N}, ::AbstractVector{N}) where {N<:Real}
remove_redundant_vertices(::VPolygon{N}; ::String="monotone_chain") where {N<:Real}
remove_redundant_vertices!(::VPolygon{N}; ::String="monotone_chain") where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N<:Real})
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractPolyhedron{N}) where {N<:Real})

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))
