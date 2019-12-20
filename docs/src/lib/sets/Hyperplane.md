```@meta
CurrentModule = LazySets
```

# [Hyperplane](@id def_Hyperplane)

```@docs
Hyperplane
dim(::Hyperplane)
ρ(::AbstractVector{N}, ::Hyperplane{N}) where {N<:Real}
σ(::AbstractVector{N}, ::Hyperplane{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Hyperplane{N}) where {N<:Real}
an_element(::Hyperplane{N}) where {N<:Real}
rand(::Type{Hyperplane})
isbounded(::Hyperplane)
isuniversal(::Hyperplane{N}, ::Bool=false) where {N<:Real}
isempty(::Hyperplane)
constrained_dimensions(::Hyperplane{N}) where {N<:Real}
constraints_list(::Hyperplane{N}) where {N<:Real}
translate(::Hyperplane{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
