```@meta
CurrentModule = LazySets
```

# [Translation](@id def_Translation)

```@docs
Translation
+(X::LazySet, v::AbstractVector)
⊕(X::LazySet, v::AbstractVector)
ρ(::AbstractVector{N}, ::Translation{N}) where {N<:Real}
σ(::AbstractVector{N}, ::Translation{N}) where {N<:Real}
an_element(::Translation)
constraints_list(::Translation{N}) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::Translation{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Translation{N}) where {N<:Real}
```
Inherited from [`AbstractAffineMap`](@ref):
* [`dim`](@ref dim(::AbstractAffineMap))
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`isbounded`](@ref isbounded(::AbstractAffineMap))
* [`vertices_list`](@ref vertices_list(::AbstractAffineMap{N}) where {N<:Real})

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
