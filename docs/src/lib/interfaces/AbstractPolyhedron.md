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

This interface requires to implement the following function:

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
constraints_list(::LazySet)
```
```@meta
CurrentModule = LazySets
```

This interface defines the following functions:

```@docs
an_element(::AbstractPolyhedron{N}) where {N}
constrained_dimensions(::AbstractPolyhedron)
isbounded(::AbstractPolyhedron{N}) where {N}
isuniversal(::AbstractPolyhedron{N}, ::Bool=false) where {N}
vertices_list(::AbstractPolyhedron)
∈(::AbstractVector, ::AbstractPolyhedron)
project(::AbstractPolyhedron{N}, ::AbstractVector{Int}) where {N}
LazySets._isbounded_stiemke
LazySets._linear_map_polyhedron
```

Some common functions implemented by several subtypes:

```@docs
addconstraint!(::AbstractPolyhedron, ::HalfSpace)
is_hyperplanar(::AbstractPolyhedron)
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
