```@meta
CurrentModule = LazySets
```

# Union

Note that the union of convex sets is generally not convex.
Hence these set types are not part of the convex-set family `ConvexSet`.

## [Binary set union (UnionSet)](@id def_UnionSet)

```@docs
UnionSet
∪(::ConvexSet, ::ConvexSet)
swap(::UnionSet)
dim(::UnionSet)
σ(::AbstractVector, ::UnionSet; algorithm="support_vector")
ρ(::AbstractVector, ::UnionSet)
an_element(::UnionSet)
∈(::AbstractVector, ::UnionSet)
isempty(::UnionSet)
isbounded(::UnionSet)
vertices_list(::UnionSet)
```

Inherited from [`ConvexSet`](@ref):
* [`singleton_list`](@ref singleton_list(::ConvexSet))

## [``n``-ary set union (UnionSetArray)](@id def_UnionSetArray)

```@docs
UnionSetArray
dim(::UnionSetArray)
array(::UnionSetArray)
σ(::AbstractVector, ::UnionSetArray; algorithm="support_vector")
ρ(::AbstractVector, ::UnionSetArray)
an_element(::UnionSetArray)
∈(::AbstractVector, ::UnionSetArray)
isempty(::UnionSetArray)
isbounded(::UnionSetArray)
vertices_list(::UnionSetArray)
```

Inherited from [`ConvexSet`](@ref):
* [`singleton_list`](@ref singleton_list(::ConvexSet))
