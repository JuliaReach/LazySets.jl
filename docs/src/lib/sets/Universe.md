```@meta
CurrentModule = LazySets
```

# [Universe](@id def_Universe)

```@docs
Universe
dim(::Universe)
ρ(::AbstractVector{N}, ::Universe{N}) where {N<:Real}
σ(::AbstractVector{N}, ::Universe{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Universe{N}) where {N<:Real}
rand(::Type{Universe})
an_element(::Universe{N}) where {N<:Real}
isempty(::Universe)
isbounded(::Universe)
isuniversal(::Universe{N}, ::Bool=false) where {N<:Real}
norm(::Universe, ::Real=Inf)
radius(::Universe, ::Real=Inf)
diameter(::Universe, ::Real=Inf)
constraints_list(::Universe{N}) where {N<:Real}
constrained_dimensions(::Universe)
translate(::Universe{N}, ::AbstractVector{N}) where {N<:Real}
```
