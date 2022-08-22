```@meta
CurrentModule = LazySets
```

# [p-norm ball (Ballp)](@id def_Ballp)

```@docs
Ballp
σ(::AbstractVector, ::Ballp)
∈(::AbstractVector, ::Ballp)
center(::Ballp)
rand(::Type{Ballp})
translate(::Ballp, ::AbstractVector)
translate!(::Ballp, ::AbstractVector)
```
Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`rectify`](@ref rectify(::ConvexSet))
* [`low`](@ref low(::ConvexSet))
* [`low`](@ref low(::ConvexSet{N}, ::Int) where {N})
* [`high`](@ref high(::ConvexSet))
* [`high`](@ref high(::ConvexSet{N}, ::Int) where {N})

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`isbounded`](@ref isbounded(::AbstractCentrallySymmetric))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetric{N}, ::Bool=false) where {N})
* [`an_element`](@ref an_element(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric, ::Int))
