```@meta
CurrentModule = LazySets
```

# [HParallelotope](@id def_HParallelotope)

```@docs
HParallelotope
directions(::HParallelotope)
offset(::HParallelotope)
dim(::HParallelotope)
base_vertex(::HParallelotope)
extremal_vertices(::HParallelotope{N, VN}) where {N, VN}
center(::HParallelotope)
genmat(::HParallelotope)
generators(::HParallelotope)
constraints_list(::HParallelotope{N, VN}) where {N, VN}
rand(::Type{HParallelotope})
isempty(P::HParallelotope)
```

Inherited from [`LazySet`](@ref):
* [`low`](@ref low(::LazySet))
* [`high`](@ref high(::LazySet))

Inherited from [`ConvexSet`](@ref):
* [`norm`](@ref norm(::ConvexSet, ::Real))
* [`radius`](@ref radius(::ConvexSet, ::Real))
* [`diameter`](@ref diameter(::ConvexSet, ::Real))
* [`singleton_list`](@ref singleton_list(::ConvexSet))
* [`low`](@ref low(::ConvexSet{N}, ::Int) where {N})
* [`high`](@ref high(::ConvexSet{N}, ::Int) where {N})
* [`rectify`](@ref rectify(::ConvexSet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`volume`](@ref volume(::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope{N}, ::Bool=false) where {N})
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope, ::Int))

Inherited from [`AbstractZonotope`](@ref):
* [`ρ`](@ref ρ(::AbstractVector, ::AbstractZonotope))
* [`σ`](@ref σ(::AbstractVector, ::AbstractZonotope))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))
* [`translate`](@ref translate(::AbstractZonotope, ::AbstractVector))
* [`vertices_list`](@ref vertices_list(::AbstractZonotope))
* [`order`](@ref order(::AbstractZonotope))
