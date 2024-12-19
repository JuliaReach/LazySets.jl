```@contents
Pages = ["AbstractSingleton.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# [Singletons (AbstractSingleton)](@id def_AbstractSingleton)

A singleton is a special hyperrectangle consisting of only one point.

```@docs
AbstractSingleton
```

This interface requires to implement the following function:

```@docs
element(::AbstractSingleton)
```

This interface defines the following functions:

```@docs
center(::AbstractSingleton)
center(::AbstractSingleton, ::Int)
element(::AbstractSingleton, ::Int)
generators(::AbstractSingleton)
genmat(::AbstractSingleton)
high(::AbstractSingleton)
high(::AbstractSingleton, ::Int)
low(::AbstractSingleton)
low(::AbstractSingleton, ::Int)
ngens(::AbstractSingleton)
radius_hyperrectangle(::AbstractSingleton)
radius_hyperrectangle(::AbstractSingleton, ::Int)
reflect(::AbstractSingleton)
vertices(::AbstractSingleton)
vertices_list(::AbstractSingleton)
∈(::AbstractVector, ::AbstractSingleton)
ρ(::AbstractVector, ::AbstractSingleton)
σ(::AbstractVector, ::AbstractSingleton)
plot_recipe(::AbstractSingleton{N}, ::Any=zero(N)) where {N}
RecipesBase.apply_recipe(::AbstractDict{Symbol,Any}, ::AbstractSingleton{N}, ::Real=zero(N)) where {N}
```

## Implementations

* [Singleton](@ref def_Singleton)
* [Origin (ZeroSet)](@ref def_ZeroSet)
