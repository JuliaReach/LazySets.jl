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
ispolyhedraltype(::Type{<:LazySet})
ispolytopic(::LazySet)
ispolytopictype(::Type{<:LazySet})
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
in(::AbstractVector, ::LazySet)
is_interior_point(::AbstractVector, ::LazySet)
linear_map(::AbstractMatrix, ::LazySet)
permute(::LazySet, ::AbstractVector{Int})
project(::LazySet, ::AbstractVector{Int})
sample(::LazySet, ::Int=1)
scale(::Real, ::LazySet)
scale!(::Real, ::LazySet)
ρ(::AbstractVector, ::LazySet)
σ(::AbstractVector, ::LazySet)
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
isapprox(::LazySet, ::LazySet)
isdisjoint(::LazySet, ::LazySet)
==(::LazySet, ::LazySet)
isequivalent(::LazySet, ::LazySet)
⊂(::LazySet, ::LazySet)
issubset(::LazySet, ::LazySet)
linear_combination(::LazySet, ::LazySet)
minkowski_difference(::LazySet, ::LazySet)
minkowski_sum(::LazySet, ::LazySet)
```
