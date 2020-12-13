```@meta
CurrentModule = LazySets
```

# [Hyperplane](@id def_Hyperplane)

```@docs
Hyperplane
dim(::Hyperplane)
ρ(::AbstractVector, ::Hyperplane)
σ(::AbstractVector, ::Hyperplane)
∈(::AbstractVector, ::Hyperplane)
an_element(::Hyperplane)
rand(::Type{Hyperplane})
isbounded(::Hyperplane)
isuniversal(::Hyperplane, ::Bool=false)
isempty(::Hyperplane)
constrained_dimensions(::Hyperplane)
constraints_list(::Hyperplane)
translate(::Hyperplane, ::AbstractVector)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))
