```@meta
CurrentModule = LazySets
```

# [Polytope in vertex representation (VPolytope)](@id def_VPolytope)

```@docs
VPolytope
dim(::VPolytope)
σ(::AbstractVector, ::VPolytope)
∈(::AbstractVector{N}, ::VPolytope{N}) where {N}
rand(::Type{VPolytope})
translate(::VPolytope, ::AbstractVector)
translate!(::VPolytope, ::AbstractVector)
vertices_list(::VPolytope)
remove_redundant_vertices(::VPolytope{N}) where {N}
constraints_list(::VPolytope)
tohrep(::VPolytope{N}) where {N}
tovrep(::VPolytope)
polyhedron(::VPolytope)
linear_map(::AbstractMatrix, ::VPolytope)
```
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`singleton_list`](@ref singleton_list(::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`linear_map`](@ref linear_map(::AbstractMatrix{NM}, ::AbstractPolyhedron{NP}) where {NM, NP})

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})
* [`volume`](@ref volume(::AbstractPolytope))
