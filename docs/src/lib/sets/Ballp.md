```@meta
CurrentModule = LazySets
```

# [p-norm ball (Ballp)](@id def_Ballp)

```@docs
Ballp
center(::Ballp)
rand(::Type{Ballp})
translate!(::Ballp, ::AbstractVector)
reflect(::Ballp)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`rectify`](@ref rectify(::LazySet))
* [`low`](@ref low(::LazySet))
* [`high`](@ref high(::LazySet))
* [`translate`](@ref translate(::LazySet, ::AbstractVector))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`isbounded`](@ref isbounded(::AbstractCentrallySymmetric))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetric{N}, ::Bool=false) where {N})
* [`an_element`](@ref an_element(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetric, ::Int))

Inherited from [`AbstractBallp`](@ref):
* [`σ`](@ref σ(::AbstractVector, ::AbstractBallp))
* [`ρ`](@ref ρ(::AbstractVector, ::AbstractBallp))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractBallp))
