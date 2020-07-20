```@meta
CurrentModule = LazySets
```

# Union

Note that the union of convex sets is generally not convex.
Hence these set types are not part of the convex-set family `LazySet`.

## [Binary set union (UnionSet)](@id def_UnionSet)

```@docs
UnionSet
∪(::LazySet, ::LazySet)
swap(::UnionSet)
dim(::UnionSet)
σ(::AbstractVector{N}, ::UnionSet{N}; algorithm="support_vector") where {N<:Real}
ρ(::AbstractVector{N}, ::UnionSet{N}) where {N<:Real}
an_element(::UnionSet{N}) where {N<:Real}
∈(::AbstractVector{N}, ::UnionSet{N}) where {N<:Real}
isempty(::UnionSet)
isbounded(::UnionSet)
vertices_list(::UnionSet)
```

Inherited from [`LazySet`](@ref):
* [`singleton_list`](@ref singleton_list(::LazySet))

## [``n``-ary set union (UnionSetArray)](@id def_UnionSetArray)

```@docs
UnionSetArray
dim(::UnionSetArray)
array(::UnionSetArray{N, S}) where {N<:Real, S<:LazySet{N}}
σ(::AbstractVector{N}, ::UnionSetArray{N}; algorithm="support_vector") where {N<:Real}
ρ(::AbstractVector{N}, ::UnionSetArray{N}) where {N<:Real}
an_element(::UnionSetArray{N}) where {N<:Real}
∈(::AbstractVector{N}, ::UnionSetArray{N}) where {N<:Real}
isempty(::UnionSetArray)
isbounded(::UnionSetArray)
vertices_list(::UnionSetArray)
```

Inherited from [`LazySet`](@ref):
* [`singleton_list`](@ref singleton_list(::LazySet))
