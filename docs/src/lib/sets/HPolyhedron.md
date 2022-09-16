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
addconstraint!(::HPoly, ::LinearConstraint)
constraints_list(::HPoly)
tohrep(::HPoly)
tovrep(::HPoly)
normalize(::HPoly{N}, p=N(2)) where {N}
isempty(::HPoly{N}, ::Bool=false) where {N}
translate(::HPoly, ::AbstractVector)
polyhedron(::HPoly)
remove_redundant_constraints(::HPoly)
remove_redundant_constraints!(::HPoly)
LazySets._isbounded_stiemke(::HPolyhedron{N}) where {N}
```

Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`singleton_list`](@ref singleton_list(::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`∈`](@ref ∈(::AbstractVector, ::AbstractPolyhedron))
* [`constrained_dimensions`](@ref constrained_dimensions(::AbstractPolyhedron)
* [`linear_map`](@ref linear_map(::AbstractMatrix{NM}, ::AbstractPolyhedron{NP}) where {NM, NP})

The following methods are specific to `HPolyhedron`.

```@docs
rand(::Type{HPolyhedron})
isbounded(::HPolyhedron)
```

Inherited from [`AbstractPolyhedron`](@ref):

* [`isuniversal`](@ref isuniversal(::AbstractPolyhedron{N}, ::Bool=false) where {N})
* [`vertices_list`](@ref vertices_list(::AbstractPolyhedron, ::Bool=false))
