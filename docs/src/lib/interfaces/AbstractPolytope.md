```@contents
Pages = ["AbstractPolytope.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Polytopes (AbstractPolytope)](@id def_AbstractPolytope)

A polytope is a bounded set with finitely many vertices (*V-representation*)
resp. facets (*H-representation*).
Note that there is a special interface combination
[Centrally symmetric polytope](@ref def_AbstractCentrallySymmetricPolytope).

```@docs
AbstractPolytope
```

This interface requires to implement the following function:

```@meta
CurrentModule = LazySets.API
```
```@docs; canonical=false
vertices_list(::LazySet)
```
```@meta
CurrentModule = LazySets
```

This interface defines the following functions:

```@docs
isbounded(::AbstractPolytope)
isempty(::AbstractPolytope)
isuniversal(::AbstractPolytope, ::Bool=false)
volume(::AbstractPolytope)
```

The following functions can be implemented for sets in vertex representation:

```@docs
remove_redundant_vertices(::AbstractPolytope)
remove_redundant_vertices!(::AbstractPolytope)
```

## Implementations

* [Polytope in constraint representation (HPolytope)](@ref def_HPolytope)
* [Polytope in vertex representation (VPolytope)](@ref def_VPolytope)
