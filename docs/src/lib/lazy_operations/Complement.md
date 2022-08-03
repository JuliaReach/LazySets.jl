```@meta
CurrentModule = LazySets
```

# [Complement](@id def_Complement)

Note that the complement of a convex set is generally not convex.
Hence this set type is not part of the convex-set family `ConvexSet`.

```@docs
Complement
dim(::Complement)
∈(::AbstractVector, ::Complement)
isempty(::Complement)
```

The concrete complement can be computed with the function `complement` (mind the lowercase,
as it is usual for functions).

```@docs
constraints_list(::Complement)
```
