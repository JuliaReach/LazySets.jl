```@meta
CurrentModule = LazySets
```

# [Half-space (HalfSpace)](@id def_HalfSpace)

```@docs
HalfSpace
LinearConstraint
dim(::HalfSpace)
ρ(::AbstractVector{N}, ::HalfSpace{N}) where {N<:Real}
σ(::AbstractVector{N}, ::HalfSpace{N}) where {N<:Real}
∈(::AbstractVector{N}, ::HalfSpace{N}) where {N<:Real}
an_element(::HalfSpace{N}) where {N<:Real}
rand(::Type{HalfSpace})
normalize(::HalfSpace{N}, p=N(2)) where {N<:Real}
isbounded(::HalfSpace)
isuniversal(::HalfSpace{N}, ::Bool=false) where {N<:Real}
isempty(::HalfSpace)
constraints_list(::HalfSpace{N}) where {N<:Real}
constraints_list(::AbstractMatrix{N}, ::AbstractVector{N}) where {N<:Real}
constrained_dimensions(::HalfSpace{N}) where {N<:Real}
translate(::HalfSpace{N}, ::AbstractVector{N}) where {N<:Real}
halfspace_left(::AbstractVector{N}, ::AbstractVector{N}) where {N<:Real}
halfspace_right(::AbstractVector{N}, ::AbstractVector{N}) where {N<:Real}
tosimplehrep(::AbstractVector{LC}) where {N<:Real, LC<:LinearConstraint{N}}
remove_redundant_constraints
remove_redundant_constraints!
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
