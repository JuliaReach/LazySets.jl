```@meta
CurrentModule = LazySets.HPolyhedronModule
```

# [Polyhedron in constraint representation (HPolyhedron)](@id def_HPolyhedron)

```@docs
HPolyhedron
```

## Operations

The following methods are shared between `HPolytope` and `HPolyhedron`.

```@docs
constraints_list(::HPoly)
dim(::HPoly)
normalize(::HPoly{N}, p::Real=N(2)) where {N}
remove_redundant_constraints(::HPoly)
remove_redundant_constraints!(::HPoly)
tohrep(::HPoly)
tovrep(::HPoly)
addconstraint!(::HPoly, ::HalfSpace)
ρ(::AbstractVector{M}, ::HPoly{N}) where {M, N}
σ(::AbstractVector{M}, ::HPoly{N}) where {M, N}
translate(::HPoly, ::AbstractVector)
```

The following methods are specific to `HPolyhedron`.

```@docs
rand(::Type{HPolyhedron})
```

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`high`](@ref high(::LazySet))
* [`isempty`](@ref isempty(::LazySet{N}, ::Bool=false) where {N})
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`reflect`](@ref reflect(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet)

Inherited from [`AbstractPolyhedron`](@ref):
* [`an_element`](@ref an_element(::AbstractPolyhedron))
* [`constrained_dimensions`](@ref constrained_dimensions(::AbstractPolyhedron))
* [`isbounded`](@ref isbounded(::AbstractPolyhedron{N}) where {N})
* [`isuniversal`](@ref isuniversal(::AbstractPolyhedron{N}, ::Bool=false) where {N})
* [`vertices_list`](@ref vertices_list(::AbstractPolyhedron, ::Bool=false))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractPolyhedron))
