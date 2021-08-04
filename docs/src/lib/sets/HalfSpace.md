```@meta
CurrentModule = LazySets
```

# [Half-space (HalfSpace)](@id def_HalfSpace)

```@docs
HalfSpace
LinearConstraint
dim(::HalfSpace)
ρ(::AbstractVector, ::HalfSpace)
σ(::AbstractVector, ::HalfSpace)
∈(::AbstractVector, ::HalfSpace)
an_element(::HalfSpace)
rand(::Type{HalfSpace})
normalize(::HalfSpace{N}, p=N(2)) where {N}
isbounded(::HalfSpace)
isuniversal(::HalfSpace, ::Bool=false)
isempty(::HalfSpace)
constraints_list(::HalfSpace)
constraints_list(::AbstractMatrix, ::AbstractVector)
constrained_dimensions(::HalfSpace)
translate(::HalfSpace, ::AbstractVector)
halfspace_left(::AbstractVector, ::AbstractVector)
halfspace_right(::AbstractVector, ::AbstractVector)
tosimplehrep(::AbstractVector{LC}) where {N, LC<:LinearConstraint{N}}
remove_redundant_constraints
remove_redundant_constraints!
complement(::HalfSpace)
project(::HalfSpace{N}, ::AbstractVector{Int}) where {N}
distance(::AbstractVector, ::HalfSpace{N}) where {N}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
