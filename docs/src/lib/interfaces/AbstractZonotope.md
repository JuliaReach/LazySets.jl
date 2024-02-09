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

This interface defines the following functions:

```@docs
ngens(::AbstractZonotope)
genmat_fallback(::AbstractZonotope{N}) where {N}
generators_fallback(::AbstractZonotope)
ρ(::AbstractVector, ::AbstractZonotope)
σ(::AbstractVector, ::AbstractZonotope)
∈(::AbstractVector, ::AbstractZonotope)
linear_map(::AbstractMatrix, ::AbstractZonotope)
translate(::AbstractZonotope, ::AbstractVector)
translate!(::AbstractZonotope, ::AbstractVector)
split(::AbstractZonotope, ::Int)
split(::AbstractZonotope, ::AbstractVector{Int}, ::AbstractVector{Int})
constraints_list(::AbstractZonotope)
constraints_list(::AbstractZonotope{N}; ::Bool=true) where {N<:AbstractFloat}
vertices_list(::AbstractZonotope)
order(::AbstractZonotope)
togrep(::AbstractZonotope)
remove_redundant_generators(::AbstractZonotope)
reduce_order(::AbstractZonotope, ::Real, ::AbstractReductionMethod=GIR05())
reflect(::AbstractZonotope)
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
