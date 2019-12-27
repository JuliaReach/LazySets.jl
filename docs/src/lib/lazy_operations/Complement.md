```@meta
CurrentModule = LazySets
```

# [Complement](@id def_Complement)

Note that the complement of a convex set is generally not convex.
Hence this set type is not part of the convex-set family `LazySet`.

```@docs
Complement
dim(::Complement)
∈(::AbstractVector{N}, ::Complement{N}) where {N<:Real}
isempty(::Complement)
```
