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

This interface defines the following functions:

```@docs
isbounded(::AbstractPolytope)
isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N}
isempty(::AbstractPolytope)
volume(::AbstractPolytope)
```

## Implementations

* [Polytope in constraint representation (HPolytope)](@ref def_HPolytope)
* [Polytope in vertex representation (VPolytope)](@ref def_VPolytope)
