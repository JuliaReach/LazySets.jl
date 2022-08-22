```@meta
CurrentModule = LazySets
```

# [Manhattan-norm ball (Ball1)](@id def_Ball1)

```@docs
Ball1
σ(::AbstractVector, ::Ball1)
∈(::AbstractVector, ::Ball1, ::Bool=false)
vertices_list(::Ball1)
center(::Ball1)
rand(::Type{Ball1})
constraints_list(::Ball1)
translate(::Ball1, ::AbstractVector)
translate!(::Ball1, ::AbstractVector)
```

Inherited from [`ConvexSet`](@ref):
* [`ρ`](@ref ρ(::AbstractVector, ::ConvexSet))
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`area`](@ref area(::ConvexSet))
* [`singleton_list`](@ref singleton_list(::ConvexSet))
* [`rectify`](@ref rectify(::ConvexSet))
* [`low`](@ref low(::ConvexSet))
* [`low`](@ref low(::ConvexSet{N}, ::Int) where {N})
* [`high`](@ref high(::ConvexSet))
* [`high`](@ref high(::ConvexSet{N}, ::Int) where {N})

Inherited from [`AbstractPolyhedron`](@ref):
* [`linear_map`](@ref linear_map(::AbstractMatrix{NM}, ::AbstractPolyhedron{NP}) where {NM, NP})

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
