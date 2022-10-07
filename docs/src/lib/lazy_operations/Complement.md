```@meta
CurrentModule = LazySets
```

# [Complement](@id def_Complement)

The concrete complement can be computed with the function `complement` (with a
lower-case "c").

Note that the complement of a convex set is generally not convex.
Hence this set type is not part of the convex-set family `ConvexSet`.

```@docs
Complement
dim(::Complement)
∈(::AbstractVector, ::Complement)
isempty(::Complement)
constraints_list(::Complement)
```
