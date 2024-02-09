```@contents
Pages = ["AbstractCentrallySymmetricPolytope.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Centrally symmetric polytopes (AbstractCentrallySymmetricPolytope)](@id def_AbstractCentrallySymmetricPolytope)

A centrally symmetric polytope is a combination of two other interfaces:
[Centrally symmetric sets](@ref def_AbstractCentrallySymmetric) and
[Polytope](@ref def_AbstractPolytope).

```@docs
AbstractCentrallySymmetricPolytope
```

This interface defines the following functions:

```@docs
dim(::AbstractCentrallySymmetricPolytope)
an_element(::AbstractCentrallySymmetricPolytope)
isempty(::AbstractCentrallySymmetricPolytope)
isuniversal(::AbstractCentrallySymmetricPolytope{N}, ::Bool=false) where {N}
center(::AbstractCentrallySymmetricPolytope, ::Int)
extrema(::AbstractCentrallySymmetricPolytope, ::Int)
extrema(::AbstractCentrallySymmetricPolytope)
```

## Implementations

* [Manhattan-norm ball (Ball1)](@ref def_Ball1)
