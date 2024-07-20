```@meta
CurrentModule = LazySets.HyperrectangleModule
```

# [Hyperrectangle](@id def_Hyperrectangle)

```@docs
Hyperrectangle
```

## Conversion

```@docs
convert(::Type{Hyperrectangle}, ::AbstractHyperrectangle)
convert(::Type{Hyperrectangle}, ::IA.IntervalBox)
convert(::Type{IA.IntervalBox}, ::AbstractHyperrectangle)
```

## Operations

```@docs
center(::Hyperrectangle)
radius_hyperrectangle(::Hyperrectangle)
radius_hyperrectangle(::Hyperrectangle, ::Int)
rand(::Type{Hyperrectangle})
permute(::Hyperrectangle, ::AbstractVector{Int})
translate(::Hyperrectangle, ::AbstractVector)
```

```@meta
CurrentModule = LazySets
```

Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))
* [`singleton_list`](@ref singleton_list(::LazySet))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope, ::Bool=false))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope))
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`order`](@ref order(::AbstractZonotope))
* [`togrep`](@ref togrep(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N})
* [`extrema`](@ref extrema(::AbstractHyperrectangle))
* [`extrema`](@ref extrema(::AbstractHyperrectangle, ::Int))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`high`](@ref high(::AbstractHyperrectangle))
* [`high`](@ref high(::AbstractHyperrectangle, ::Int))
* [`isflat`](@ref isflat(::Hyperrectangle))
* [`low`](@ref low(::AbstractHyperrectangle))
* [`low`](@ref low(::AbstractHyperrectangle, ::Int))
* [`ngens`](@ref ngens(::AbstractHyperrectangle))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`rectify`](@ref rectify(::AbstractHyperrectangle))
* [`reflect`](@ref reflect(::AbstractHyperrectangle))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle))
* [`∈`](@ref ∈(::AbstractVector, ::AbstractHyperrectangle))
* [`ρ`](@ref ρ(::AbstractVector, ::AbstractHyperrectangle))
* [`σ`](@ref σ(::AbstractVector, ::AbstractHyperrectangle))
