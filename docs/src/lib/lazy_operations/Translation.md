```@meta
CurrentModule = LazySets
```

# [Translation](@id def_Translation)

```@docs
Translation
+(X::LazySet, v::AbstractVector)
⊕(X::LazySet, v::AbstractVector)
ρ(::AbstractVector, ::Translation)
σ(::AbstractVector, ::Translation)
an_element(::Translation)
constraints_list(::Translation)
linear_map(::AbstractMatrix, ::Translation)
∈(::AbstractVector, ::Translation)
center(::Translation)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`reflect`](@ref reflect(::LazySet))

Inherited from [`AbstractAffineMap`](@ref):
* [`dim`](@ref dim(::AbstractAffineMap))
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`isbounded`](@ref isbounded(::AbstractAffineMap))
* [`vertices_list`](@ref vertices_list(::AbstractAffineMap))
