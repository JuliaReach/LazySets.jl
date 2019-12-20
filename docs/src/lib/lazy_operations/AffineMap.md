```@meta
CurrentModule = LazySets
```

# [Affine map (AffineMap)](@id def_AffineMap)

```@docs
AffineMap
dim(::AffineMap)
σ(::AbstractVector{N}, ::AffineMap{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::AffineMap{N}) where {N<:Real}
an_element(::AffineMap)
isempty(::AffineMap)
isbounded(::AffineMap)
∈(::AbstractVector{N}, ::AffineMap{N}) where {N<:Real}
vertices_list(::AffineMap{N}) where {N<:Real}
constraints_list(::AffineMap{N}) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::AffineMap{N}) where {N<:Real}
```
