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

Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`singleton_list`](@ref singleton_list(::ConvexSet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N})
* [`volume`](@ref volume(::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`∈`](@ref ∈(::AbstractVector, ::AbstractZonotope))
* [`translate`](@ref translate(::AbstractZonotope, ::AbstractVector))
* [`constraints_list`](@ref constraints_list(::AbstractZonotope{N}; ::Bool=true) where {N<:AbstractFloat})
* [`order`](@ref order(::AbstractZonotope))
