```@meta
CurrentModule = LazySets
```

# [Polytope in constraint representation (HPolytope)](@id def_HPolytope)

[Convex polytopes](https://en.wikipedia.org/wiki/Polytope) are bounded polyhedra.
The type `HPolytope` represents polytopes.
While identical to [`HPolyhedron`](@ref) in implementation, `HPolytope`
instances are assumed to be bounded.

```@docs
HPolytope
```

Some functionality is shared with [`HPolyhedron`](@ref).
Below follows the additional functionality specific to `HPolytope`.

```@docs
rand(::Type{HPolytope})
vertices_list(::HPolytope{N}) where {N<:Real}
isbounded(::HPolytope, ::Bool=true)
```

Inherited from [`AbstractPolytope`](@ref):
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N<:Real})
