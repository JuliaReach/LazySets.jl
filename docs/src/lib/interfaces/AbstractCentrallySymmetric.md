```@contents
Pages = ["AbstractCentrallySymmetric.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Centrally symmetric sets (AbstractCentrallySymmetric)](@id def_AbstractCentrallySymmetric)

Centrally symmetric sets such as balls of different norms are characterized by a
center.
Note that there is a special interface combination
[Centrally symmetric polytope](@ref def_AbstractCentrallySymmetricPolytope).

```@docs
AbstractCentrallySymmetric
```

This interface defines the following functions:

```@docs
dim(::AbstractCentrallySymmetric)
isbounded(::AbstractCentrallySymmetric)
isuniversal(::AbstractCentrallySymmetric{N}, ::Bool=false) where {N}
an_element(::AbstractCentrallySymmetric)
isempty(::AbstractCentrallySymmetric)
center(::AbstractCentrallySymmetric, ::Int)
extrema(::AbstractCentrallySymmetric, ::Int)
extrema(::AbstractCentrallySymmetric)
```

## Implementations

* [Euclidean-norm ball (Ball2)](@ref def_Ball2)
* [Ellipsoid](@ref def_Ellipsoid)
* [p-norm ball (Ballp)](@ref def_Ballp)
