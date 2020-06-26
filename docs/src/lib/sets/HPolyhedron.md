```@meta
CurrentModule = LazySets
```

# [Polyhedron in constraint representation (HPolyhedron)](@id def_HPolyhedron)

```@docs
HPolyhedron
```

The following methods are shared between `HPolytope` and `HPolyhedron`.

```@docs
dim(::HPoly{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::HPoly{N}) where {N<:Real}
σ(::AbstractVector{N}, ::HPoly{N}) where {N<:Real}
addconstraint!(::HPoly{N}, ::LinearConstraint{N}) where {N<:Real}
constraints_list(::HPoly{N}) where {N<:Real}
tohrep(::HPoly{N}) where {N<:Real}
tovrep(::HPoly{N}) where {N<:Real}
normalize(::HPoly{N}, p=N(2)) where {N<:Real}
isempty(::HPoly{N}, ::Bool=false) where {N<:Real}
translate(::HPoly{N}, ::AbstractVector{N}) where {N<:Real}
polyhedron(::HPoly{N}) where {N<:Real}
remove_redundant_constraints(::HPoly{N}) where {N<:Real}
remove_redundant_constraints!(::HPoly{N}) where {N<:Real}
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolyhedron`](@ref):
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractPolyhedron{N}) where {N<:Real})
* [`constrained_dimensions`](@ref constrained_dimensions(::AbstractPolyhedron)
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractPolyhedron{N}) where {N<:Real})

The following methods are specific to `HPolyhedron`.

```@docs
rand(::Type{HPolyhedron})
isbounded(::HPolyhedron)
vertices_list(::HPolyhedron{N}) where {N<:Real}
singleton_list(::HPolyhedron{N}) where {N<:Real}
```

Inherited from [`AbstractPolyhedron`](@ref):
* [`isuniversal`](@ref isuniversal(::AbstractPolyhedron{N}, ::Bool=false) where {N<:Real})
