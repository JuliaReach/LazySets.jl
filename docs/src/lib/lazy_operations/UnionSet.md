```@meta
CurrentModule = LazySets
```

# Union

## [Binary set union (UnionSet)](@id def_UnionSet)

```@docs
UnionSet
∪(::LazySet, ::LazySet)
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

## [``n``-ary set union (UnionSetArray)](@id def_UnionSetArray)

```@docs
UnionSetArray
UnionSet!
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
