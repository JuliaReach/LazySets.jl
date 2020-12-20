```@meta
CurrentModule = LazySets
```

# [Euclidean-norm ball (Ball2)](@id def_Ball2)

```@docs
Ball2
ρ(::AbstractVector, ::Ball2)
σ(::AbstractVector, ::Ball2)
∈(::AbstractVector, ::Ball2)
center(::Ball2)
rand(::Type{Ball2})
sample(::Ball2{N}, ::Int) where {N}
translate(::Ball2, ::AbstractVector)
chebyshev_center(::Ball2)
volume(::Ball2)
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
