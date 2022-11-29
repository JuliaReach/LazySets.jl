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
constraints_list(::HParallelotope)
rand(::Type{HParallelotope})
isempty(::HParallelotope)
volume(::HParallelotope)
```

Inherited from [`LazySet`](@ref):
* [`low`](@ref low(::LazySet))
* [`high`](@ref high(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))
* [`rectify`](@ref rectify(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))

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
