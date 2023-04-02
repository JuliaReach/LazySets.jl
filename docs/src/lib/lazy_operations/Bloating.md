```@meta
CurrentModule = LazySets
```

# [Bloating](@id def_Bloating)

```@docs
Bloating
dim(::Bloating)
σ(::AbstractVector, ::Bloating)
ρ(::AbstractVector, ::Bloating)
isbounded(::Bloating)
isempty(::Bloating)
an_element(::Bloating)
constraints_list(::Bloating)
center(::Bloating)
is_polyhedral(::Bloating)
```
Inherited from [`LazySet`](@ref):
* [`high`](@ref high(::LazySet))
* [`low`](@ref low(::LazySet))
* [`reflect`](@ref reflect(::LazySet))
