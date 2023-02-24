```@meta
CurrentModule = LazySets
```

# [Polyhedron in constraint representation (HPolyhedron)](@id def_HPolyhedron)

```@docs
HPolyhedron
```

The following methods are shared between `HPolytope` and `HPolyhedron`.

```@docs
dim(::HPoly)
ρ(::AbstractVector{M}, ::HPoly{N}) where {M, N}
σ(::AbstractVector{M}, ::HPoly{N}) where {M, N}
addconstraint!(::HPoly, ::HalfSpace)
constraints_list(::HPoly)
tohrep(::HPoly)
tovrep(::HPoly)
normalize(::HPoly{N}, p=N(2)) where {N}
translate(::HPoly, ::AbstractVector)
remove_redundant_constraints(::HPoly)
remove_redundant_constraints!(::HPoly)
```
Inherited from [`LazySet`](@ref):
* [`high`](@ref high(::LazySet))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`isempty`](@ref isempty(::LazySet{N}, ::Bool=false) where {N})

Inherited from [`AbstractPolyhedron`](@ref):
* [`∈`](@ref ∈(::AbstractVector, ::AbstractPolyhedron))
* [`an_element`](@ref an_element(::AbstractPolyhedron))
* [`constrained_dimensions`](@ref constrained_dimensions(::AbstractPolyhedron))
* [`linear_map`](@ref linear_map(::AbstractMatrix{NM}, ::AbstractPolyhedron{NP}) where {NM, NP})

The following methods are specific to `HPolyhedron`.

```@docs
rand(::Type{HPolyhedron})
```

Inherited from [`AbstractPolyhedron`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolyhedron{N}) where {N})
* [`isuniversal`](@ref isuniversal(::AbstractPolyhedron{N}, ::Bool=false) where {N})
* [`vertices_list`](@ref vertices_list(::AbstractPolyhedron, ::Bool=false))
