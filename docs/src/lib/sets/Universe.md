```@meta
CurrentModule = LazySets
```

# [Universe](@id def_Universe)

```@docs
Universe
dim(::Universe)
ρ(::AbstractVector, ::Universe)
σ(::AbstractVector, ::Universe)
∈(::AbstractVector, ::Universe)
rand(::Type{Universe})
an_element(::Universe{N}) where {N}
isempty(::Universe)
isbounded(::Universe)
isuniversal(::Universe{N}, ::Bool=false) where {N}
norm(::Universe, ::Real=Inf)
radius(::Universe, ::Real=Inf)
diameter(::Universe, ::Real=Inf)
constraints(::Universe{N}) where {N}
constraints_list(::Universe{N}) where {N}
constrained_dimensions(::Universe)
translate(::Universe, ::AbstractVector)
translate!(::Universe, ::AbstractVector)
permute(::Universe, ::AbstractVector{Int})
complement(::Universe{N}) where {N}
polyhedron(::Universe)
```
