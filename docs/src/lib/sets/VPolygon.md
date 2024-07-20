```@meta
CurrentModule = LazySets.VPolygonModule
```

# [Polygon in vertex representation (VPolygon)](@id def_VPolygon)

```@docs
VPolygon
```

## Conversion

```@docs
convert(::Type{VPolygon}, ::LazySet)
```

## Operations

```@docs
area(::VPolygon)
an_element(::VPolygon)
constraints_list(::VPolygon)
rand(::Type{VPolygon})
remove_redundant_vertices(::VPolygon; ::String="monotone_chain")
remove_redundant_vertices!(::VPolygon; ::String="monotone_chain")
tohrep(::VPolygon{N}, ::Type{HPOLYGON}=HPolygon) where {N, HPOLYGON<:AbstractHPolygon}
tovrep(::VPolygon)
vertices_list(::VPolygon)
∈(::AbstractVector, ::VPolygon)
linear_map(::AbstractMatrix, ::VPolygon)
permute(::VPolygon, ::AbstractVector{Int})
σ(::AbstractVector, ::VPolygon)
translate(::VPolygon, ::AbstractVector)
translate!(::VPolygon, ::AbstractVector)
```

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`reflect`](@ref reflect(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))
* [`volume`](@ref volume(::AbstractPolygon))
