```@meta
CurrentModule = LazySets
```

# [Polytope in constraint representation (HPolytope)](@id def_HPolytope)

[Convex polytopes](https://en.wikipedia.org/wiki/Polytope) are bounded polyhedra.
The type `HPolytope` represents polytopes.
While identical to [`HPolyhedron`](@ref) in its representation, `HPolytope`
instances are assumed to be bounded.

```@docs
HPolytope
```

Most functionality is shared with [`HPolyhedron`](@ref).
Additional functionality specific to `HPolytope` is listed below.

```@docs
rand(::Type{HPolytope})
vertices_list(::HPolytope)
isbounded(::HPolytope, ::Bool=true)
```

Inherited from [`AbstractPolytope`](@ref):
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})
* [`volume`](@ref volume(::AbstractPolytope))
