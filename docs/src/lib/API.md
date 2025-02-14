```@contents
Pages = ["API.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets.API
```

```@docs
API
```

# Set interface

```@docs
LazySet
```

# Unary set functions

```@docs
an_element(::LazySet)
area(::LazySet)
center(::LazySet, ::Int)
center(::LazySet)
complement(::LazySet)
concretize(::LazySet)
constraints_list(::LazySet)
constraints(::LazySet)
convex_hull(::LazySet)
diameter(::LazySet, ::Real=Inf)
dim(::LazySet)
eltype(::Type{<:LazySet})
eltype(::LazySet)
extrema(::LazySet, ::Int)
extrema(::LazySet)
high(::LazySet, ::Int)
high(::LazySet)
isbounded(::LazySet)
isboundedtype(::Type{<:LazySet})
isconvextype(::Type{<:LazySet})
isempty(::LazySet, ::Bool=false)
isoperation(::LazySet)
isoperationtype(::Type{<:LazySet})
ispolyhedral(::LazySet)
isuniversal(::LazySet, ::Bool=false)
low(::LazySet, ::Int)
low(::LazySet)
norm(::LazySet, ::Real=Inf)
radius(::LazySet, ::Real=Inf)
rand(::Type{<:LazySet})
rectify(::LazySet)
reflect(::LazySet)
vertices_list(::LazySet)
vertices(::LazySet)
volume(::LazySet)
```

# Mixed set functions

```@docs
affine_map(::AbstractMatrix, ::LazySet, ::AbstractVector)
distance(::AbstractVector, ::LazySet)
exponential_map(::AbstractMatrix, ::LazySet)
∈(::AbstractVector, ::LazySet)
is_interior_point(::AbstractVector, ::LazySet)
linear_map(::AbstractMatrix, ::LazySet)
permute(::LazySet, ::AbstractVector{Int})
project(::LazySet, ::AbstractVector{Int})
sample(::LazySet, ::Int=1)
scale(::Real, ::LazySet)
scale!(::Real, ::LazySet)
ρ(::AbstractVector, ::LazySet)
support_function
σ(::AbstractVector, ::LazySet)
support_vector
translate(::LazySet, ::AbstractVector)
translate!(::LazySet, ::AbstractVector)
```

# Binary set functions

```@docs
cartesian_product(::LazySet, ::LazySet)
convex_hull(::LazySet, ::LazySet)
difference(::LazySet, ::LazySet)
distance(::LazySet, ::LazySet)
exact_sum(::LazySet, ::LazySet)
intersection(::LazySet, ::LazySet)
≈(::LazySet, ::LazySet)
isdisjoint(::LazySet, ::LazySet)
is_intersection_empty
==(::LazySet, ::LazySet)
isequivalent(::LazySet, ::LazySet)
⊂(::LazySet, ::LazySet)
⊆(::LazySet, ::LazySet)
linear_combination(::LazySet, ::LazySet)
minkowski_difference(::LazySet, ::LazySet)
pontryagin_difference
minkowski_sum(::LazySet, ::LazySet)
```
