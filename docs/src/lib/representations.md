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
σ(::AbstractVector{AbstractFloat}, ::Ball2{AbstractFloat})
∈(::AbstractVector{AbstractFloat}, ::Ball2{AbstractFloat})
center(::Ball2)
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
center(::BallInf)
radius(::BallInf, ::Real)
radius_hyperrectangle(::BallInf)
radius_hyperrectangle(::BallInf, ::Int)
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
σ(::AbstractVector{Real}, ::Ball1{Real})
∈(::AbstractVector{Real}, ::Ball1{Real})
vertices_list(::Ball1)
center(::Ball1)
constraints_list(::Ball1)
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
σ(::AbstractVector{AbstractFloat}, ::Ballp{AbstractFloat})
∈(::AbstractVector{AbstractFloat}, ::Ballp{AbstractFloat})
center(::Ballp)
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
σ(::AbstractVector{AbstractFloat}, ::Ellipsoid{AbstractFloat})
∈(::AbstractVector{AbstractFloat}, ::Ellipsoid{AbstractFloat})
center(::Ellipsoid)
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
σ(::AbstractVector{Real}, ::EmptySet{Real})
∈(::AbstractVector{Real}, ::EmptySet{Real})
an_element(::EmptySet)
isempty(::EmptySet)
norm(::EmptySet, ::Real)
radius(::EmptySet, ::Real)
diameter(::EmptySet, ::Real)
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
σ(::AbstractVector{Real}, ::HalfSpace{Real})
∈(::AbstractVector{Real}, ::HalfSpace{Real})
an_element(::HalfSpace{N}) where {N<:Real}
isempty(::HalfSpace)
constraints_list(::HalfSpace{N}) where {N<:Real}
constrained_dimensions(::HalfSpace{N}) where {N<:Real}
LazySets.halfspace_left(::AbstractVector{Real}, ::AbstractVector{Real})
LazySets.halfspace_right(::AbstractVector{Real}, ::AbstractVector{Real})
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

## Hyperplane

```@docs
Hyperplane
dim(::Hyperplane)
σ(::AbstractVector{Real}, ::Hyperplane{Real})
∈(::AbstractVector{Real}, ::Hyperplane{Real})
an_element(::Hyperplane{N}) where {N<:Real}
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
center(::Hyperrectangle)
radius_hyperrectangle(::Hyperrectangle)
radius_hyperrectangle(::Hyperrectangle, ::Int)
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
σ(::AbstractVector{Real}, ::Interval{Real})
∈(::AbstractVector, ::Interval)
∈(::Real, ::Interval)
an_element(::Interval)
vertices_list(::Interval)
center(::Interval)
low(::Interval)
high(::Interval)
+(::Interval, ::Interval)
-(::Interval, ::Interval)
*(::Interval, ::Interval)
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
σ(::AbstractVector{Real}, ::Line{Real})
∈(::AbstractVector{Real}, ::Line{Real})
an_element(::Line{N}) where {N<:Real}
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
σ(::AbstractVector{Real}, ::LineSegment{Real})
∈(::AbstractVector{Real}, ::LineSegment{Real})
LazySets.halfspace_left(::LineSegment)
LazySets.halfspace_right(::LineSegment)
LazySets.constraints_list(::LineSegment)
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
σ(::AbstractVector{Real}, ::HPolygon{Real})
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
σ(::AbstractVector{Real}, ::HPolygonOpt{Real})
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
σ(::AbstractVector{Real}, ::VPolygon{Real})
∈(::AbstractVector{Real}, ::VPolygon{Real})
an_element(::VPolygon{N}) where {N<:Real}
vertices_list(::VPolygon)
tohrep(::VPolygon)
tovrep(::VPolygon)
constraints_list(::VPolygon)
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
<=(::AbstractVector{AbstractFloat}, ::AbstractVector{AbstractFloat})
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

The following methods are shared between polyhedra and polytopes. 

```@docs
dim(::HPoly{Real})
ρ(::AbstractVector{Real}, ::HPoly{Real})
σ(::AbstractVector{Real}, ::HPoly{Real})
∈(::AbstractVector{Real}, ::HPoly{Real})
addconstraint!(::HPoly{Real}, ::LinearConstraint{Real})
constraints_list(::HPoly{Real})
tosimplehrep(::HPoly{Real})
tohrep(::HPoly{Real})
isempty(::HPoly{N}) where {N<:Real}
convex_hull(P1::HPoly{Real}, P2::HPoly{Real})
cartesian_product(P1::HPoly{Real}, P2::HPoly{Real})
tovrep(::HPoly{Real})
vertices_list(::HPolyhedron{Real})
singleton_list(::HPolyhedron{N}) where {N<:Real}
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

### Vertex representation

```@docs
VPolytope
dim(::VPolytope)
σ(::AbstractVector{Real}, ::VPolytope{Real})
vertices_list(::VPolytope)
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
element(::Singleton)
element(::Singleton, ::Int)
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
element(::ZeroSet)
element(::ZeroSet, ::Int)
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
σ(::AbstractVector{Real}, ::Zonotope{Real})
∈(::AbstractVector{Real}, ::Zonotope{Real})
vertices_list(::Zonotope)
center(::Zonotope)
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
