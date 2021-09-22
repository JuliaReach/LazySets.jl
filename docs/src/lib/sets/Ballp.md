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
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`isbounded`](@ref isbounded(::AbstractCentrallySymmetric))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetric{N}, ::Bool=false) where {N})
* [`an_element`](@ref an_element(::AbstractCentrallySymmetric))
