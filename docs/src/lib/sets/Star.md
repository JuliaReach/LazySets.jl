```@meta
CurrentModule = LazySets
```

# [Star](@id def_Star)

```@docs
Star
center(::Star)
predicate(::Star)
basis(::Star)
dim(::Star)
σ(::AbstractVector, ::Star)
ρ(::AbstractVector, ::Star)
an_element(::Star)
isempty(::Star)
isbounded(::Star)
∈(::AbstractVector, ::Star)
vertices_list(::Star)
constraints_list(::Star)
linear_map(::AbstractMatrix, ::Star)
affine_map(::AbstractMatrix, ::Star, v::AbstractVector)
rand(::Type{Star})
```
Inherited from [`LazySet`](@ref):
* [`high`](@ref high(::LazySet))
* [`low`](@ref low(::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`singleton_list`](@ref singleton_list(::ConvexSet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`isuniversal`](@ref isuniversal(::AbstractPolyhedron{N}, ::Bool=false) where {N})
