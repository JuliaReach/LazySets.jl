```@contents
Pages = ["AbstractZonotope.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Zonotopes (AbstractZonotope)](@id def_AbstractZonotope)

A zonotope is a specific centrally symmetric polytope characterized by a
center and a collection of generators.

```@docs
AbstractZonotope
```

This interface requires to implement the following functions:

```@docs
generators(::AbstractZonotope)
genmat(::AbstractZonotope)
```

This interface defines the following functions:

```@docs
constraints_list(::AbstractZonotope)
constraints_list(::AbstractZonotope{N}; ::Bool=true) where {N<:AbstractFloat}
generators_fallback(::AbstractZonotope)
genmat_fallback(::AbstractZonotope{N}) where {N}
ngens(::AbstractZonotope)
order(::AbstractZonotope)
reflect(::AbstractZonotope)
remove_redundant_generators(::AbstractZonotope)
togrep(::AbstractZonotope)
vertices_list(::AbstractZonotope)
∈(::AbstractVector, ::AbstractZonotope)
linear_map(::AbstractMatrix, ::AbstractZonotope)
reduce_order(::AbstractZonotope, ::Real, ::AbstractReductionMethod=GIR05())
split(::AbstractZonotope, ::Int)
split(::AbstractZonotope, ::AbstractVector{Int}, ::AbstractVector{Int})
ρ(::AbstractVector, ::AbstractZonotope)
σ(::AbstractVector, ::AbstractZonotope)
```

## Order-reduction methods

```@docs
LazySets.AbstractReductionMethod
LazySets.ASB10
LazySets.COMB03
LazySets.GIR05
```

## Implementations

* [Zonotope](@ref def_Zonotope)
* [Line segment (LineSegment)](@ref def_LineSegment)
