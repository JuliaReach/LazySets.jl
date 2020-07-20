```@meta
CurrentModule = LazySets
```

# [Empty set (EmptySet)](@id def_EmptySet)

```@docs
EmptySet
∅
dim(::EmptySet)
σ(::AbstractVector{N}, ::EmptySet{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::EmptySet{N}) where {N<:Real}
∈(::AbstractVector{N}, ::EmptySet{N}) where {N<:Real}
an_element(::EmptySet)
rand(::Type{EmptySet})
isbounded(::EmptySet)
isempty(::EmptySet)
isuniversal(::EmptySet{N}, ::Bool=false) where {N<:Real}
norm(::EmptySet, ::Real=Inf)
radius(::EmptySet, ::Real=Inf)
diameter(::EmptySet, ::Real=Inf)
linear_map(::AbstractMatrix{N}, ::EmptySet{N}) where {N}
translate(::EmptySet{N}, ::AbstractVector{N}) where {N<:Real}
plot_recipe(::EmptySet{N}, ::N=zero(N)) where {N<:Real}
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::EmptySet{N}, ::N=zero(N)) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))
