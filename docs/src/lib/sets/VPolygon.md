```@meta
CurrentModule = LazySets
```

# [Polygon in vertex representation (VPolygon)](@id def_VPolygon)

```@docs
VPolygon
σ(::AbstractVector, ::VPolygon)
∈(::AbstractVector, ::VPolygon)
an_element(::VPolygon)
rand(::Type{VPolygon})
vertices_list(::VPolygon)
tohrep(::VPolygon{N}, ::Type{HPOLYGON}=HPolygon) where {N, HPOLYGON<:AbstractHPolygon}
tovrep(::VPolygon)
constraints_list(::VPolygon)
translate(::VPolygon, ::AbstractVector)
translate!(::VPolygon, ::AbstractVector)
remove_redundant_vertices(::VPolygon; ::String="monotone_chain")
remove_redundant_vertices!(::VPolygon; ::String="monotone_chain")
permute(::VPolygon, ::AbstractVector{Int})
linear_map(M::AbstractMatrix, P::VPolygon)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`reflect`](@ref reflect(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})
* [`volume`](@ref volume(::AbstractPolytope))

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))
