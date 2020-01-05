```@meta
CurrentModule = LazySets
```

# [Manhattan-norm ball (Ball1)](@id def_Ball1)

```@docs
Ball1
σ(::AbstractVector{N}, ::Ball1{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Ball1{N}, ::Bool=false) where {N<:Real}
vertices_list(::Ball1{N, VN}) where {N<:Real, VN<:AbstractVector{N}}
center(::Ball1{N}) where {N<:Real}
rand(::Type{Ball1})
constraints_list(::Ball1{N}) where {N<:Real}
translate(::Ball1{N}, ::AbstractVector{N}) where {N<:Real}
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isuniversal`](@ref isuniversal(::AbstractPolytope{N}, ::Bool=false) where {N<:Real})
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope{N}) where {N<:Real})
