```@meta
CurrentModule = LazySets.StarModule
```

# [Star](@id def_Star)

```@docs
Star
```

## Conversion

```@docs
convert(::Type{Star}, ::AbstractPolyhedron{N}) where {N}
```

## Operations

```@docs
an_element(::Star)
basis(::Star)
center(::Star)
constraints_list(::Star)
dim(::Star)
isempty(::Star)
isbounded(::Star)
predicate(::Star)
rand(::Type{Star})
vertices_list(::Star)
affine_map(::AbstractMatrix, ::Star, v::AbstractVector)
∈(::AbstractVector, ::Star)
linear_map(::AbstractMatrix, ::Star)
ρ(::AbstractVector, ::Star)
σ(::AbstractVector, ::Star)
```

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`high`](@ref high(::LazySet))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`reflect`](@ref reflect(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolyhedron`](@ref):
* [`isuniversal`](@ref isuniversal(::AbstractPolyhedron{N}, ::Bool=false) where {N})
