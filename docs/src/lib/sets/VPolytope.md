```@meta
CurrentModule = LazySets.VPolytopeModule
```

# [Polytope in vertex representation (VPolytope)](@id def_VPolytope)

```@docs
VPolytope
```

## Conversion

```@docs
convert(::Type{VPolytope}, ::LazySet)
```

## Operations

```@docs
constraints_list(::VPolytope)
dim(::VPolytope)
polyhedron(::VPolytope)
rand(::Type{VPolytope})
remove_redundant_vertices(::VPolytope{N}) where {N}
reflect(::VPolytope)
tohrep(::VPolytope{N}) where {N}
tovrep(::VPolytope)
vertices_list(::VPolytope)
linear_map(::AbstractMatrix, ::VPolytope)
∈(::AbstractVector{N}, ::VPolytope{N}) where {N}
ρ(::AbstractVector, ::VPolytope)
σ(::AbstractVector, ::VPolytope)
translate(::VPolytope, ::AbstractVector)
translate!(::VPolytope, ::AbstractVector)
cartesian_product(::VPolytope, ::VPolytope)
```

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`high`](@ref high(::LazySet))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`an_element`](@ref an_element(::AbstractPolyhedron))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})
* [`volume`](@ref volume(::AbstractPolytope))
