```@meta
CurrentModule = LazySets
```

# [Translation](@id def_Translation)

```@docs
Translation
+(X::ConvexSet, v::AbstractVector)
⊕(X::ConvexSet, v::AbstractVector)
ρ(::AbstractVector, ::Translation)
σ(::AbstractVector, ::Translation)
an_element(::Translation)
constraints_list(::Translation)
linear_map(::AbstractMatrix, ::Translation)
∈(::AbstractVector, ::Translation)
center(::Translation)
```
Inherited from [`AbstractAffineMap`](@ref):
* [`dim`](@ref dim(::AbstractAffineMap))
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`isbounded`](@ref isbounded(::AbstractAffineMap))
* [`vertices_list`](@ref vertices_list(::AbstractAffineMap))

Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`singleton_list`](@ref singleton_list(::ConvexSet))
