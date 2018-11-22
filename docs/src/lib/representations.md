# Common Set Representations

This section of the manual describes the basic set representation types.

```@contents
Pages = ["representations.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
DocTestSetup = quote
    using LazySets
    using Compat.SparseArrays, Compat.LinearAlgebra
end
```

## Balls

### Euclidean norm ball

```@docs
Ball2
σ(::AbstractVector{N}, ::Ball2{N}) where {N<:AbstractFloat}
∈(::AbstractVector{N}, ::Ball2{N}) where {N<:AbstractFloat}
center(::Ball2{N}) where {N<:AbstractFloat}
rand(::Type{Ball2})
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* `an_element`

### Infinity norm ball

```@docs
BallInf
center(::BallInf{N}) where {N<:Real}
radius(::BallInf, ::Real=Inf)
radius_hyperrectangle(::BallInf{N}) where {N<:Real}
radius_hyperrectangle(::BallInf{N}, ::Int) where {N<:Real}
rand(::Type{BallInf})
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* `an_element`

Inherited from [`AbstractHyperrectangle`](@ref):
* [`σ`](@ref σ(::AbstractVector{Real}, ::AbstractHyperrectangle{Real}))
* [`∈`](@ref ∈(::AbstractVector{Real}, ::AbstractHyperrectangle{Real}))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle{Real}))
* [`high`](@ref high(::AbstractHyperrectangle{Real}))
* [`low`](@ref low(::AbstractHyperrectangle{Real}))

### Manhattan norm ball

```@docs
Ball1
σ(::AbstractVector{N}, ::Ball1{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Ball1{N}) where {N<:Real}
vertices_list(::Ball1{N}) where {N<:Real}
center(::Ball1{N}) where {N<:Real}
rand(::Type{Ball1})
constraints_list(::Ball1{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* `an_element`

### p-norm ball

```@docs
Ballp
σ(::AbstractVector{N}, ::Ballp{N}) where {N<:AbstractFloat}
∈(::AbstractVector{N}, ::Ballp{N}) where {N<:AbstractFloat}
center(::Ballp{N}) where {N<:AbstractFloat}
rand(::Type{Ballp})
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* `an_element`

## Ellipsoid

```@docs
Ellipsoid
σ(::AbstractVector{N}, ::Ellipsoid{N}) where {N<:AbstractFloat}
∈(::AbstractVector{N}, ::Ellipsoid{N}) where {N<:AbstractFloat}
rand(::Type{Ellipsoid})
center(::Ellipsoid{N}) where {N<:AbstractFloat}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* `an_element`

## Empty set

```@docs
EmptySet
∅
dim(::EmptySet)
σ(::AbstractVector{N}, ::EmptySet{N}) where {N<:Real}
∈(::AbstractVector{N}, ::EmptySet{N}) where {N<:Real}
an_element(::EmptySet)
rand(::Type{EmptySet})
isempty(::EmptySet)
norm(::EmptySet, ::Real=Inf)
radius(::EmptySet, ::Real=Inf)
diameter(::EmptySet, ::Real=Inf)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

## Half-Space

```@docs
HalfSpace
LinearConstraint
dim(::HalfSpace)
σ(::AbstractVector{N}, ::HalfSpace{N}) where {N<:Real}
∈(::AbstractVector{N}, ::HalfSpace{N}) where {N<:Real}
an_element(::HalfSpace{N}) where {N<:Real}
rand(::Type{HalfSpace})
isempty(::HalfSpace)
constraints_list(::HalfSpace{N}) where {N<:Real}
constrained_dimensions(::HalfSpace{N}) where {N<:Real}
halfspace_left(::AbstractVector{N}, ::AbstractVector{N}) where {N<:Real}
halfspace_right(::AbstractVector{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

## Hyperplane

```@docs
Hyperplane
dim(::Hyperplane)
σ(::AbstractVector{N}, ::Hyperplane{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Hyperplane{N}) where {N<:Real}
an_element(::Hyperplane{N}) where {N<:Real}
rand(::Type{Hyperplane})
isempty(::Hyperplane)
constrained_dimensions(::Hyperplane{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

## Hyperrectangle

```@docs
Hyperrectangle
rand(::Type{Hyperrectangle})
center(::Hyperrectangle{N}) where {N<:Real}
radius_hyperrectangle(::Hyperrectangle{N}) where {N<:Real}
radius_hyperrectangle(::Hyperrectangle{N}, ::Int) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* `an_element`

Inherited from [`AbstractHyperrectangle`](@ref):
* [`σ`](@ref σ(::AbstractVector{Real}, ::AbstractHyperrectangle{Real}))
* [`∈`](@ref ∈(::AbstractVector{Real}, ::AbstractHyperrectangle{Real}))
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle{Real}))
* [`high`](@ref high(::AbstractHyperrectangle{Real}))
* [`low`](@ref low(::AbstractHyperrectangle{Real}))

## Interval

```@docs
Interval
dim(::Interval)
σ(::AbstractVector{N}, ::Interval{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Interval{N}) where {N<:Real}
∈(::N, ::Interval{N}) where {N<:Real}
an_element(::Interval{N}) where {N<:Real}
vertices_list(::Interval{N}) where {N<:Real}
center(::Interval{N}) where {N<:Real}
low(::Interval)
high(::Interval)
radius_hyperrectangle(::Interval{N}) where {N<:Real}
radius_hyperrectangle(::Interval{N}, ::Int) where {N<:Real}
+(::Interval, ::Interval)
-(::Interval, ::Interval)
*(::Interval, ::Interval)
rand(::Type{Interval})
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))

## Line

```@docs
Line
dim(::Line)
σ(::AbstractVector{N}, ::Line{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Line{N}) where {N<:Real}
an_element(::Line{N}) where {N<:Real}
rand(::Type{Line})
isempty(::Line)
constrained_dimensions(::Line{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

## Line segment

```@docs
LineSegment
dim(::LineSegment)
σ(::AbstractVector{N}, ::LineSegment{N}) where {N<:Real}
∈(::AbstractVector{N}, ::LineSegment{N}) where {N<:Real}
center(::LineSegment{N}) where {N<:Real}
an_element(::LineSegment{N}) where {N<:Real}
rand(::Type{LineSegment})
halfspace_left(::LineSegment)
halfspace_right(::LineSegment)
vertices_list(::LineSegment{N}) where {N<:Real}
constraints_list(::LineSegment{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

## Polygons

### Constraint representation

```@docs
HPolygon
σ(::AbstractVector{N}, ::HPolygon{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))

Inherited from [`AbstractHPolygon`](@ref):
* `an_element`
* [`∈`](@ref ∈(::AbstractVector{Real}, ::AbstractHPolygon{Real}))
* [`vertices_list`](@ref vertices_list(::AbstractHPolygon{Real}))
* [`tohrep`](@ref tohrep(::AbstractHPolygon{Real}))
* [`tovrep`](@ref tovrep(::AbstractHPolygon{Real}))
* [`addconstraint!`](@ref addconstraint!(::AbstractHPolygon{Real}, ::LinearConstraint{Real}))
* [`constraints_list`](@ref constraints_list(::AbstractHPolygon{Real}))

### Optimized constraint representation

```@docs
HPolygonOpt
σ(::AbstractVector{N}, ::HPolygonOpt{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))

Inherited from [`AbstractHPolygon`](@ref):
* `an_element`
* [`∈`](@ref ∈(::AbstractVector{Real}, ::AbstractHPolygon{Real}))
* [`vertices_list`](@ref vertices_list(::AbstractHPolygon{Real}))
* [`tohrep`](@ref tohrep(::AbstractHPolygon{Real}))
* [`tovrep`](@ref tovrep(::AbstractHPolygon{Real}))
* [`addconstraint!`](@ref addconstraint!(::AbstractHPolygon{Real}, ::LinearConstraint{Real}))
* [`constraints_list`](@ref constraints_list(::AbstractHPolygon{Real}))

### Vertex representation

```@docs
VPolygon
σ(::AbstractVector{N}, ::VPolygon{N}) where {N<:Real}
∈(::AbstractVector{N}, ::VPolygon{N}) where {N<:Real}
an_element(::VPolygon{N}) where {N<:Real}
rand(::Type{VPolygon})
vertices_list(::VPolygon{N}) where {N<:Real}
tohrep(::VPolygon{N}, ::Type{HPOLYGON}=HPolygon) where {N<:Real, HPOLYGON<:AbstractHPolygon}
tovrep(::VPolygon{N}) where {N<:Real}
constraints_list(::VPolygon{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))

### Sorting directions

```@docs
LazySets.jump2pi
<=(::AbstractVector{N}, ::AbstractVector{N}) where {N<:AbstractFloat}
<=(::AbstractVector{N}, ::AbstractVector{N}) where {N<:Real}
LazySets.quadrant(::AbstractVector{Real})
```
## Polyhedra and Polytopes


### Constraint representation

[Convex polytopes](https://en.wikipedia.org/wiki/Polytope) are bounded polyhedra.
The types `HPolytope` and `HPolyhedron` are used to represent polytopes and
general polyhedra respectively, the difference being that for `HPolytope` there
is a running assumption about the boundedness of the set.

```@docs
HPolytope
HPolyhedron
```

The following methods are shared between `HPolytope` and `HPolyhedron`.

```@docs
dim(::HPoly{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::HPoly{N}) where {N<:Real}
σ(::AbstractVector{N}, ::HPoly{N}) where {N<:Real}
∈(::AbstractVector{N}, ::HPoly{N}) where {N<:Real}
addconstraint!(::HPoly{N}, ::LinearConstraint{N}) where {N<:Real}
constraints_list(::HPoly{N}) where {N<:Real}
copy(P::PT) where {N, PT<:HPoly{N}} where {N<:Real}
tosimplehrep(::HPoly{N}) where {N<:Real}
tohrep(::HPoly{N}) where {N<:Real}
isempty(::HPoly{N}) where {N<:Real}
convex_hull(::HPoly{N}, ::HPoly{N}) where {N<:Real}
cartesian_product(::HPoly{N}, ::HPoly{N}) where {N<:Real}
tovrep(::HPoly{N}) where {N<:Real}
polyhedron(::HPoly{N}) where {N<:Real}
remove_redundant_constraints
remove_redundant_constraints!
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

#### Polytopes in constraint representation

The following methods are specific for `HPolytope`.

```@docs
rand(::Type{HPolytope})
vertices_list(::HPolytope{N}) where {N<:Real}
```

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
The following methods are specific for polytopes.

#### Polyhedra

The following methods are specific for `HPolyhedron`.

```@docs
rand(::Type{HPolyhedron})
vertices_list(::HPolyhedron{N}) where {N<:Real}
singleton_list(::HPolyhedron{N}) where {N<:Real}
```

### Vertex representation

```@docs
VPolytope
dim(::VPolytope)
σ(::AbstractVector{N}, ::VPolytope{N}) where {N<:Real}
rand(::Type{VPolytope})
vertices_list(::VPolytope{N}) where {N<:Real}
constraints_list(::VPolytope{N}) where {N<:Real}
tohrep(::VPolytope{N}) where {N<:Real}
tovrep(::VPolytope)
cartesian_product(::VPolytope{N}, ::VPolytope{N}) where N
polyhedron(::VPolytope{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

## Singleton

```@docs
Singleton
rand(::Type{Singleton})
element(::Singleton{N}) where {N<:Real}
element(::Singleton{N}, ::Int) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`high`](@ref high(::AbstractHyperrectangle{Real}))
* [`low`](@ref low(::AbstractHyperrectangle{Real}))

Inherited from [`AbstractSingleton`](@ref):
* `σ`
* `∈`
* `an_element`
* [`center`](@ref center(::AbstractSingleton{Real}))
* `vertices_list`
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{Real}))
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{Real}, ::Int))
* `linear_map`

## Zero set

```@docs
ZeroSet
dim(::ZeroSet)
σ(::AbstractVector{N}, ::ZeroSet{N}) where {N<:Real}
∈(::AbstractVector{N}, ::ZeroSet{N}) where {N<:Real}
rand(::Type{ZeroSet})
element(::ZeroSet{N}) where {N<:Real}
element(::ZeroSet{N}, ::Int) where {N<:Real}
linear_map(::AbstractMatrix, ::ZeroSet{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`high`](@ref high(::AbstractHyperrectangle{Real}))
* [`low`](@ref low(::AbstractHyperrectangle{Real}))

Inherited from [`AbstractSingleton`](@ref):
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{Real}))
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{Real}, ::Int))
* `vertices_list`
* [`center`](@ref center(::AbstractSingleton{Real}))
* `an_element`

## Zonotope

```@docs
Zonotope
σ(::AbstractVector{N}, ::Zonotope{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Zonotope{N}) where {N<:Real}
rand(::Type{Zonotope})
vertices_list(::Zonotope{N}) where {N<:Real}
constraints_list(::Zonotope{N}) where {N<:Real}
center(::Zonotope{N}) where {N<:Real}
order(::Zonotope)
minkowski_sum(::Zonotope, ::Zonotope)
linear_map(::AbstractMatrix, ::Zonotope)
scale(::Real, ::Zonotope)
ngens(::Zonotope)
reduce_order(::Zonotope, r)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* `an_element`
