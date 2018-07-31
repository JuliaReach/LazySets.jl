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

Inherited from [`AbstractPointSymmetric`](@ref):
* [`dim`](@ref dim(::AbstractPointSymmetric))
* `an_element`

### Infinity norm ball

```@docs
BallInf
dim(::BallInf)
σ(::AbstractVector{Real}, ::BallInf{Real})
∈(::AbstractVector{Real}, ::BallInf{Real})
an_element(::BallInf)
vertices_list(::BallInf)
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

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))

### Manhattan norm ball

```@docs
Ball1
dim(::Ball1)
σ(::AbstractVector{Real}, ::Ball1{Real})
∈(::AbstractVector{Real}, ::Ball1{Real})
an_element(::Ball1)
vertices_list(::Ball1)
center(::Ball1)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

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

Inherited from [`AbstractPointSymmetric`](@ref):
* [`dim`](@ref dim(::AbstractPointSymmetric))
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

Inherited from [`AbstractPointSymmetric`](@ref):
* [`dim`](@ref dim(::AbstractPointSymmetric))
* `an_element`

## Empty set

```@docs
EmptySet
∅
dim(::EmptySet)
σ(::AbstractVector{Real}, ::EmptySet{Real})
∈(::AbstractVector{Real}, ::EmptySet{Real})
an_element(::EmptySet)
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
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

## Hyperrectangle

```@docs
Hyperrectangle
Hyperrectangle(;kwargs...)
dim(::Hyperrectangle)
σ(::AbstractVector{Real}, ::Hyperrectangle{Real})
∈(::AbstractVector{Real}, ::Hyperrectangle{Real})
an_element(::Hyperrectangle)
vertices_list(::Hyperrectangle)
center(::Hyperrectangle)
radius_hyperrectangle(::Hyperrectangle)
radius_hyperrectangle(::Hyperrectangle, ::Int)
high(::Hyperrectangle)
low(::Hyperrectangle)
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))

## Interval

```@docs
Interval
dim(::Interval)
σ(::AbstractVector{Real}, ::Interval{Real, IntervalArithmetic.AbstractInterval{Real}})
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
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

## Polygons

### Constraint representation

```@docs
HPolygon
dim(::HPolygon)
σ(::AbstractVector{Real}, ::HPolygon{Real})
∈(::AbstractVector{Real}, ::HPolygon{Real})
an_element(::HPolygon)
vertices_list(::HPolygon)
tohrep(::HPolygon)
tovrep(::HPolygon)
addconstraint!(::HPolygon{Real}, ::LinearConstraint{Real})
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

### Optimized constraint representation

```@docs
HPolygonOpt
dim(::HPolygonOpt)
σ(::AbstractVector{Real}, ::HPolygonOpt{Real})
∈(::AbstractVector{Real}, ::HPolygonOpt{Real})
an_element(::HPolygonOpt)
vertices_list(::HPolygonOpt)
tohrep(::HPolygonOpt)
tovrep(::HPolygonOpt)
addconstraint!(::HPolygonOpt{Real}, ::LinearConstraint{Real})
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
VPolygon
dim(::VPolygon)
σ(::AbstractVector{Real}, ::VPolygon{Real})
∈(::AbstractVector{Real}, ::VPolygon{Real})
an_element(::VPolygon{N}) where {N<:Real}
vertices_list(::VPolygon)
tohrep(::VPolygon)
tovrep(::VPolygon)
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

### Sorting directions

```@docs
LazySets.jump2pi
<=(::AbstractVector{AbstractFloat}, ::AbstractVector{AbstractFloat})
LazySets.quadrant(::AbstractVector{Real})
```

## Polytopes

### Constraint representation

```@docs
HPolytope
dim(::HPolytope)
σ(::AbstractVector{Real}, ::HPolytope{Real})
∈(::AbstractVector{Real}, ::HPolytope{Real})
addconstraint!(::HPolytope{Real}, ::LinearConstraint{Real})
constraints_list(::HPolytope)
tosimplehrep(::HPolytope)
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
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

## Singleton

```@docs
Singleton
dim(::Singleton)
σ(::AbstractVector{Real}, ::Singleton{Real})
∈(::AbstractVector{Real}, ::Singleton{Real})
an_element(::Singleton)
vertices_list(::Singleton)
center(::Singleton)
radius_hyperrectangle(::Singleton)
radius_hyperrectangle(::Singleton, ::Int)
element(::Singleton)
element(::Singleton, ::Int)
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))
* [`linear_map`](@ref linear_map(::AbstractMatrix, ::AbstractPolytope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))

## Zero set

```@docs
ZeroSet
dim(::ZeroSet)
σ(::AbstractVector{N}, ::ZeroSet) where {N<:Real}
∈(::AbstractVector{N}, ::ZeroSet{N}) where {N<:Real}
an_element(::ZeroSet)
vertices_list(::ZeroSet)
center(::ZeroSet)
radius_hyperrectangle(::ZeroSet)
radius_hyperrectangle(::ZeroSet, ::Int)
element(::ZeroSet)
element(::ZeroSet, ::Int)
linear_map(::AbstractMatrix, ::ZeroSet{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))

## Zonotope

```@docs
Zonotope
dim(::Zonotope)
σ(::AbstractVector{Real}, ::Zonotope{Real})
∈(::AbstractVector{Real}, ::Zonotope{Real})
an_element(::Zonotope)
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
