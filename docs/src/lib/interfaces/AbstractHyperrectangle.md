```@contents
Pages = ["AbstractHyperrectangle.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Hyperrectangles (AbstractHyperrectangle)](@id def_AbstractHyperrectangle)

A hyperrectangle is a special centrally symmetric polytope with axis-aligned
facets.

```@docs
AbstractHyperrectangle
```

This interface defines the following functions:

```@docs
□(c, r)
norm(::AbstractHyperrectangle, ::Real=Inf)
radius(::AbstractHyperrectangle, ::Real=Inf)
σ(::AbstractVector, ::AbstractHyperrectangle)
ρ(::AbstractVector, ::AbstractHyperrectangle)
∈(::AbstractVector, ::AbstractHyperrectangle)
vertices_list(::AbstractHyperrectangle)
constraints_list(::AbstractHyperrectangle{N}) where {N}
high(::AbstractHyperrectangle)
high(::AbstractHyperrectangle, ::Int)
low(::AbstractHyperrectangle)
low(::AbstractHyperrectangle, ::Int)
extrema(::AbstractHyperrectangle)
extrema(::AbstractHyperrectangle, ::Int)
isflat(::AbstractHyperrectangle)
split(::AbstractHyperrectangle{N}, ::AbstractVector{Int}) where {N}
generators(::AbstractHyperrectangle)
genmat(::AbstractHyperrectangle)
ngens(::AbstractHyperrectangle{N}) where {N}
rectify(::AbstractHyperrectangle)
volume(::AbstractHyperrectangle)
distance(::AbstractVector, ::AbstractHyperrectangle{N}; ::Real=N(2)) where {N}
reflect(::AbstractHyperrectangle)
```

## Implementations

Concrete set representations:

* [Hyperrectangle](@ref def_Hyperrectangle)
* [Infinity-norm ball (BallInf)](@ref def_BallInf)
* [Interval](@ref def_Interval)

Lazy set operations:

* [Symmetric interval hull (SymmetricIntervalHull)](@ref def_SymmetricIntervalHull)
