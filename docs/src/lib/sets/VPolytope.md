```@meta
CurrentModule = LazySets
```

# [Polytope in vertex representation (VPolytope)](@id def_VPolytope)

```@docs
VPolytope
dim(::VPolytope)
σ(::AbstractVector{N}, ::VPolytope{N}) where {N<:Real}
∈(::AbstractVector{N}, ::VPolytope{N}) where {N<:Real}
rand(::Type{VPolytope})
translate(::VPolytope{N}, ::AbstractVector{N}) where {N<:Real}
vertices_list(::VPolytope{N}) where {N<:Real}
remove_redundant_vertices(::VPolytope{N}) where {N<:Real}
constraints_list(::VPolytope{N}) where {N<:Real}
tohrep(::VPolytope{N}) where {N<:Real}
tovrep(::VPolytope)
polyhedron(::VPolytope{N}) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::VPolytope{N}) where {N<:Real}
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
