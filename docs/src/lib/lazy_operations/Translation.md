```@meta
CurrentModule = LazySets
```

# [Translation](@id def_Translation)

```@docs
Translation
+(X::LazySet, v::AbstractVector)
⊕(X::LazySet, v::AbstractVector)
dim(::Translation)
ρ(::AbstractVector{N}, ::Translation{N}) where {N<:Real}
σ(::AbstractVector{N}, ::Translation{N}) where {N<:Real}
an_element(::Translation)
isempty(::Translation)
constraints_list(::Translation{N}) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::Translation{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Translation{N}) where {N<:Real}
```
