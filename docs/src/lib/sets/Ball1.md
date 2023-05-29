```@meta
CurrentModule = LazySets
```

# [Manhattan-norm ball (Ball1)](@id def_Ball1)

```@docs
Ball1
σ(::AbstractVector, ::Ball1)
ρ(::AbstractVector, ::Ball1)
∈(::AbstractVector, ::Ball1, ::Bool=false)
vertices_list(::Ball1)
center(::Ball1)
rand(::Type{Ball1})
constraints_list(::Ball1)
translate(::Ball1, ::AbstractVector)
translate!(::Ball1, ::AbstractVector)
reflect(::Ball1)
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`area`](@ref area(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`rectify`](@ref rectify(::LazySet))
* [`low`](@ref low(::LazySet))
* [`high`](@ref high(::LazySet))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::LazySet)

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`volume`](@ref volume(::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope{N}, ::Bool=false) where {N})
* [`center`](@ref center(::AbstractCentrallySymmetricPolytope, ::Int))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope, ::Int))
