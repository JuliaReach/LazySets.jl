```@meta
CurrentModule = LazySets
```

# [Linear map (LinearMap)](@id def_LinearMap)

```@docs
LinearMap
*(::Union{AbstractMatrix, UniformScaling, AbstractVector, Real}, ::ConvexSet)
dim(::LinearMap)
ρ(::AbstractVector, ::LinearMap)
σ(::AbstractVector, ::LinearMap)
∈(::AbstractVector, ::LinearMap)
an_element(::LinearMap)
vertices_list(::LinearMap)
constraints_list(::LinearMap)
linear_map(::AbstractMatrix, ::LinearMap)
project(S::ConvexSet{N}, ::AbstractVector{Int}, ::Type{LM}, ::Int=dim(S)) where {N, LM<:LinearMap}
```
Inherited from [`AbstractAffineMap`](@ref):
* [`isempty`](@ref isempty(::AbstractAffineMap))
* [`isbounded`](@ref isbounded(::AbstractAffineMap))

Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`singleton_list`](@ref singleton_list(::ConvexSet))

The lazy projection of a set can be conveniently constructed using `Projection`.

```@docs
Projection
```
