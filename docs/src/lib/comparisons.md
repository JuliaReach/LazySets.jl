# Comparisons

This section of the manual lists the comparison functions in floating point between scalars and between vectors.

```@contents
Pages = ["comparisons.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

## Tolerance type

```@docs
Tolerance
```

## Approximate inequality

```@docs
_leq(x::N, y::N; kwargs...) where {N<:Real}
_leq(x::N, y::M; kwargs...) where {N<:Real, M<:Real}
_geq(x::Real, y::Real; kwargs...)
_leq(x::N, y::N; rtol::Real=Base.rtoldefault(N), ztol::Real=ABSZTOL(N), atol::Real=zero(N)) where {N<:AbstractFloat}
```

## Approximate equality

```@docs
_isapprox(x::N, y::N; rtol::Real=Base.rtoldefault(N), ztol::Real=ABSZTOL(N), atol::Real=zero(N)) where {N<:AbstractFloat}
isapproxzero(x::Real; kwargs...)
isapproxzero(x::N; ztol::Real=ABSZTOL(N)) where {N<:AbstractFloat}
```
