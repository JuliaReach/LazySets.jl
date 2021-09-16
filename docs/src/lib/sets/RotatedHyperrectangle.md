```@meta
CurrentModule = LazySets
```

# [Rotated hyperrectangle (RotatedHyperrectangle)](@id def_RotatedHyperrectangle)

```@docs
RotatedHyperrectangle
dim(::RotatedHyperrectangle)
ρ(::AbstractVector, ::RotatedHyperrectangle)
σ(::AbstractVector, ::RotatedHyperrectangle)
center(::RotatedHyperrectangle)
generators(::RotatedHyperrectangle)
genmat(::RotatedHyperrectangle)
linear_map(::AbstractMatrix, ::RotatedHyperrectangle)
vertices_list(::RotatedHyperrectangle)
constraints_list(::RotatedHyperrectangle)
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`∈`](@ref ∈(::AbstractVector, ::AbstractZonotope))
* [`translate`](@ref translate(::AbstractZonotope, ::AbstractVector))
* [`constraints_list`](@ref constraints_list(::AbstractZonotope{N}; ::Bool=true) where {N<:AbstractFloat})
* [`order`](@ref order(::AbstractZonotope))
