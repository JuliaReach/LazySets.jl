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
reflect(::VPolytope)
```
Inherited from [`LazySet`](@ref):
* [`high`](@ref high(::LazySet))
* [`low`](@ref low(::LazySet))
* [`ρ`](@ref ρ(::AbstractVector, ::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`an_element`](@ref an_element(::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})
* [`volume`](@ref volume(::AbstractPolytope))
