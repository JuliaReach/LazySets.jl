```@meta
CurrentModule = LazySets
```

# [Bloating](@id def_Bloating)

```@docs
Bloating
dim(::Bloating)
σ(::AbstractVector{N}, ::Bloating{N}) where {N<:AbstractFloat}
ρ(::AbstractVector{N}, ::Bloating{N}) where {N<:AbstractFloat}
isbounded(::Bloating)
isempty(::Bloating)
an_element(::Bloating)
```

Inherited from [`LazySet`](@ref):
* [`singleton_list`](@ref singleton_list(::LazySet))
