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
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`ρ`](@ref ρ(::AbstractVector, ::AbstractZonotope))
* [`σ`](@ref σ(::AbstractVector, ::AbstractZonotope))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))
* [`translate`](@ref translate(::AbstractZonotope, ::AbstractVector))
* [`constraints_list`](@ref constraints_list(::AbstractZonotope))
* [`constraints_list`](@ref constraints_list(::AbstractZonotope{N}; ::Bool=true) where {N<:AbstractFloat})
* [`vertices_list`](@ref vertices_list(::AbstractZonotope))
* [`order`](@ref order(::AbstractZonotope))
