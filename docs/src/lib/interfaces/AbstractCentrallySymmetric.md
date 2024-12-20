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

This interface requires to implement the following functions:

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
center(::LazySet)
```
```@meta
CurrentModule = LazySets
```

This interface defines the following functions:

```@docs
an_element(::AbstractCentrallySymmetric)
dim(::AbstractCentrallySymmetric)
center(::AbstractCentrallySymmetric, ::Int)
extrema(::AbstractCentrallySymmetric)
extrema(::AbstractCentrallySymmetric, ::Int)
isbounded(::AbstractCentrallySymmetric)
isempty(::AbstractCentrallySymmetric)
isuniversal(::AbstractCentrallySymmetric, ::Bool=false)
```

## Implementations

* [Euclidean-norm ball (Ball2)](@ref def_Ball2)
* [Ellipsoid](@ref def_Ellipsoid)
* [p-norm ball (Ballp)](@ref def_Ballp)
