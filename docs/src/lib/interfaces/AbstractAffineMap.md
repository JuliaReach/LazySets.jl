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
an_element(::AbstractAffineMap)
center(::AbstractAffineMap)
constraints_list(::AbstractAffineMap)
dim(::AbstractAffineMap)
isbounded(::AbstractAffineMap)
isempty(::AbstractAffineMap)
∈(::AbstractVector, ::AbstractAffineMap)
linear_map(::AbstractMatrix, ::AbstractAffineMap)
ρ(::AbstractVector, ::AbstractAffineMap)
σ(::AbstractVector, ::AbstractAffineMap)
vertices_list(::AbstractAffineMap)
```

## Implementations

* [Affine map (AffineMap)](@ref def_AffineMap)
* [Exponential map (ExponentialMap)](@ref def_ExponentialMap)
* [Linear map (LinearMap)](@ref def_LinearMap)
* [Reset map (ResetMap)](@ref def_ResetMap)
* [Translation](@ref def_Translation)
