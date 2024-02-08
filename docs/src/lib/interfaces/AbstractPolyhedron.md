```@contents
Pages = ["AbstractPolyhedron.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Polyhedra (AbstractPolyhedron)](@id def_AbstractPolyhedron)

A polyhedron has finitely many facets (*H-representation*) and is not
necessarily bounded.

```@docs
AbstractPolyhedron
```

This interface defines the following functions:

```@docs
∈(::AbstractVector, ::AbstractPolyhedron)
isuniversal(::AbstractPolyhedron{N}, ::Bool=false) where {N}
constrained_dimensions(::AbstractPolyhedron)
an_element(::AbstractPolyhedron{N}) where {N}
isbounded(::AbstractPolyhedron{N}) where {N}
vertices_list(::AbstractPolyhedron)
project(::AbstractPolyhedron{N}, ::AbstractVector{Int}) where {N}
LazySets._linear_map_polyhedron
LazySets._isbounded_stiemke
```

Some common functions to work with linear constraints:

```@docs
constraints_list(::AbstractMatrix, ::AbstractVector)
tosimplehrep(::AbstractVector{LC}) where {N, LC<:HalfSpace{N}}
remove_redundant_constraints(::AbstractVector{S}) where {S<:HalfSpace}
remove_redundant_constraints!(::AbstractVector{S}) where {S<:HalfSpace}
```

Plotting (bounded) polyhedra is available, too:

```@docs
plot_recipe(::AbstractPolyhedron{N}, ::Any=zero(N)) where {N}
```

## Implementations

* [Half-space (HalfSpace)](@ref def_HalfSpace)
* [Polyhedron in constraint representation (HPolyhedron)](@ref def_HPolyhedron)
* [Hyperplane](@ref def_Hyperplane)
* [Line2D](@ref def_Line2D)
* [Line](@ref def_Line)
* [Universe](@ref def_Universe)
* [Star](@ref def_Star)
