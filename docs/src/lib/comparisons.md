# Comparisons

This section of the manual lists the comparison functions in floating point between scalars and between vectors.

```@contents
Pages = ["comparisons.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
DocTestSetup = quote
    using LazySets
end
```

```@docs
_leq(x::N, y::N) where {N<:Real}
_leq(x::N, y::M) where {N<:Real, M<:Real}
_geq(x::N, y::N) where {N<:Real}
_geq(x::N, y::M) where {N<:Real, M<:Real}
isapproxzero(x::N; kwargs...) where {N<:Real}
isapproxzero(x::N; ztol::Real=ABSZTOL(N)) where {N<:AbstractFloat}
_isapprox(x::N, y::N; rtol::Real=Base.rtoldefault(N), ztol::Real=ABSZTOL(N), atol::Real=zero(N)) where {N<:AbstractFloat}
_leq(x::N, y::N; rtol::Real=Base.rtoldefault(N), ztol::Real=ABSZTOL(N), atol::Real=zero(N)) where {N<:AbstractFloat}
```
