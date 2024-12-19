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

This interface requires to implement the following function:

```@docs
radius_hyperrectangle(::AbstractHyperrectangle)
```

This interface defines the following functions:

```@docs
□(c, r)
constraints_list(::AbstractHyperrectangle)
extrema(::AbstractHyperrectangle)
extrema(::AbstractHyperrectangle, ::Int)
generators(::AbstractHyperrectangle)
genmat(::AbstractHyperrectangle)
high(::AbstractHyperrectangle)
high(::AbstractHyperrectangle, ::Int)
isflat(::AbstractHyperrectangle)
low(::AbstractHyperrectangle)
low(::AbstractHyperrectangle, ::Int)
ngens(::AbstractHyperrectangle{N}) where {N}
norm(::AbstractHyperrectangle, ::Real=Inf)
radius(::AbstractHyperrectangle, ::Real=Inf)
radius_hyperrectangle(::AbstractHyperrectangle, ::Int)
rectify(::AbstractHyperrectangle)
reflect(::AbstractHyperrectangle)
vertices_list(::AbstractHyperrectangle)
volume(::AbstractHyperrectangle)
distance(::AbstractVector, ::AbstractHyperrectangle{N}; ::Real=N(2)) where {N}
∈(::AbstractVector, ::AbstractHyperrectangle)
split(::AbstractHyperrectangle{N}, ::AbstractVector{Int}) where {N}
ρ(::AbstractVector, ::AbstractHyperrectangle)
σ(::AbstractVector, ::AbstractHyperrectangle)
```

## Implementations

Concrete set representations:

* [Hyperrectangle](@ref def_Hyperrectangle)
* [Infinity-norm ball (BallInf)](@ref def_BallInf)
* [Interval](@ref def_Interval)

Lazy set operations:

* [Symmetric interval hull (SymmetricIntervalHull)](@ref def_SymmetricIntervalHull)
