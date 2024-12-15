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
an_element(::AbstractCentrallySymmetricPolytope)
dim(::AbstractCentrallySymmetricPolytope)
center(::AbstractCentrallySymmetricPolytope, ::Int)
extrema(::AbstractCentrallySymmetricPolytope)
extrema(::AbstractCentrallySymmetricPolytope, ::Int)
isempty(::AbstractCentrallySymmetricPolytope)
isuniversal(::AbstractCentrallySymmetricPolytope, ::Bool=false)
```

## Implementations

* [Manhattan-norm ball (Ball1)](@ref def_Ball1)
