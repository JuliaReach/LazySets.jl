```@meta
CurrentModule = LazySets.HParallelotopeModule
```

# [HParallelotope](@id def_HParallelotope)

```@docs
HParallelotope
```

## Conversion

```@docs
convert(::Type{HParallelotope}, Z::AbstractZonotope{N}) where {N}
```

## Operations

```@docs
base_vertex(::HParallelotope)
center(::HParallelotope)
constraints_list(::HParallelotope)
dim(::HParallelotope)
directions(::HParallelotope)
extremal_vertices(::HParallelotope{N, VN}) where {N, VN}
generators(::HParallelotope)
genmat(::HParallelotope)
offset(::HParallelotope)
rand(::Type{HParallelotope})
volume(::HParallelotope)
```

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`high`](@ref high(::LazySet))
* [`low`](@ref low(::LazySet))
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`rectify`](@ref rectify(::LazySet))
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope))
* [`extrema`](@ref extrema(::AbstractCentrallySymmetricPolytope, ::Int))
* [`isuniversal`](@ref isuniversal(::AbstractCentrallySymmetricPolytope, ::Bool=false))

Inherited from [`AbstractZonotope`](@ref):
* [`isempty`](@ref isempty(::AbstractZonotope))
* [`order`](@ref order(::AbstractZonotope))
* [`reflect`](@ref reflect(::AbstractZonotope))
* [`vertices_list`](@ref vertices_list(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractZonotope))
* [`ρ`](@ref ρ(::AbstractVector, ::AbstractZonotope))
* [`σ`](@ref σ(::AbstractVector, ::AbstractZonotope))
* [`translate`](@ref translate(::AbstractZonotope, ::AbstractVector))
