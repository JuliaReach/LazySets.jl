```@meta
CurrentModule = LazySets
```

# [Empty set (EmptySet)](@id def_EmptySet)

```@docs
EmptySet
∅
dim(::EmptySet)
σ(::AbstractVector, ::EmptySet)
ρ(::AbstractVector, ::EmptySet)
∈(::AbstractVector, ::EmptySet)
an_element(::EmptySet)
rand(::Type{EmptySet})
isbounded(::EmptySet)
isempty(::EmptySet)
isuniversal(::EmptySet{N}, ::Bool=false) where {N}
vertices(::EmptySet{N}) where {N}
vertices_list(::EmptySet{N}) where {N}
norm(::EmptySet, ::Real=Inf)
radius(::EmptySet, ::Real=Inf)
diameter(::EmptySet, ::Real=Inf)
linear_map(::AbstractMatrix{N}, ::EmptySet{N}) where {N}
translate(::EmptySet, ::AbstractVector)
plot_recipe(::EmptySet{N}, ::Any=zero(N)) where {N}
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::EmptySet{N}, ::N=zero(N)) where {N}
area(::EmptySet{N}) where {N}
```
Inherited from [`LazySet`](@ref):
* [`singleton_list`](@ref singleton_list(::LazySet))
