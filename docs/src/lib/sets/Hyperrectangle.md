```@meta
CurrentModule = LazySets
```

# [Hyperrectangle](@id def_Hyperrectangle)

```@docs
Hyperrectangle
rand(::Type{Hyperrectangle})
center(::Hyperrectangle)
radius_hyperrectangle(::Hyperrectangle)
radius_hyperrectangle(::Hyperrectangle, ::Int)
translate(::Hyperrectangle, ::AbstractVector)
permute(::Hyperrectangle, ::AbstractVector{Int})
```
Inherited from [`ConvexSet`](@ref):
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`singleton_list`](@ref singleton_list(::ConvexSet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope, ::Bool=false))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`σ`](@ref σ(::AbstractVector, ::AbstractHyperrectangle))
* [`ρ`](@ref ρ(::AbstractVector, ::AbstractHyperrectangle))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractHyperrectangle))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle))
* [`low`](@ref low(::AbstractHyperrectangle))
* [`low`](@ref low(::AbstractHyperrectangle, ::Int))
* [`high`](@ref high(::AbstractHyperrectangle))
* [`high`](@ref high(::AbstractHyperrectangle, ::Int))
* [`extrema`](@ref extrema(::AbstractHyperrectangle))
* [`extrema`](@ref extrema(::AbstractHyperrectangle, ::Int))
* [`isflat`](@ref isflat(::Hyperrectangle))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`ngens`](@ref ngens(::AbstractHyperrectangle))
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N})
* [`rectify`](@ref rectify(::AbstractHyperrectangle))
