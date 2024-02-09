```@contents
Pages = ["AbstractAffineMap.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Affine maps (AbstractAffineMap)](@id def_AbstractAffineMap)

An affine map consists of a linear map and a translation.

```@docs
AbstractAffineMap
```

This interface defines the following functions:

```@docs
dim(::AbstractAffineMap)
σ(::AbstractVector, ::AbstractAffineMap)
ρ(::AbstractVector, ::AbstractAffineMap)
an_element(::AbstractAffineMap)
isempty(::AbstractAffineMap)
isbounded(::AbstractAffineMap)
∈(::AbstractVector, ::AbstractAffineMap)
center(::AbstractAffineMap)
vertices_list(::AbstractAffineMap)
constraints_list(::AbstractAffineMap)
linear_map(::AbstractMatrix, ::AbstractAffineMap)
```

## Implementations

* [Affine map (AffineMap)](@ref def_AffineMap)
* [Exponential map (ExponentialMap)](@ref def_ExponentialMap)
* [Linear map (LinearMap)](@ref def_LinearMap)
* [Reset map (ResetMap)](@ref def_ResetMap)
* [Translation](@ref def_Translation)
