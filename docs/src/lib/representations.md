# Common Set Representations

This section of the manual describes the basic set representation types.

```@contents
Pages = ["representations.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

## Abstract support function and support vector

```@docs
LazySets
LazySets.LazySet
ρ
```

## Balls

### Euclidean norm ball

```@docs
Ball2
dim(B::Ball2)
σ(d::AbstractVector{Float64}, B::Ball2)
```

### Infinity norm ball

```@docs
BallInf
dim(B::BallInf)
σ(d::AbstractVector{Float64}, B::BallInf)
vertices_list(B::BallInf)
norm(B::BallInf, p::Real=Inf)
radius(B::BallInf, p::Real=Inf)
diameter(B::BallInf, p::Real=Inf)
```

## Polygons

### Constraint representation

```@docs
HPolygon
addconstraint!(P::HPolygon, c::LinearConstraint)
dim(P::HPolygon)
σ(d::AbstractVector{Float64}, P::HPolygon)

is_contained(x::AbstractVector{Float64}, P::HPolygon)
tovrep(P::HPolygon)
vertices_list(P::HPolygon)
```

### Optimized constraint representation

```@docs
HPolygonOpt
addconstraint!(P::HPolygonOpt, c::LinearConstraint)
dim(P::HPolygonOpt)
σ(d::AbstractVector{Float64}, P::HPolygonOpt)

is_contained(x::AbstractVector{Float64}, P::HPolygonOpt)
tovrep(P::HPolygonOpt)
vertices_list(P::HPolygonOpt)
```

### Vertex representation

```@docs
VPolygon
dim(P::VPolygon)
vertices_list(P::VPolygon)
singleton_list(P::VPolygon)
```

## Lines and linear constraints

```@docs
intersection
LinearConstraint
Line
```

## Hyperrectangles

```@docs
Hyperrectangle
dim(H::Hyperrectangle)
σ(d::AbstractVector{Float64}, H::Hyperrectangle)
vertices_list(H::Hyperrectangle)
norm(H::Hyperrectangle, p::Real=Inf)
radius(H::Hyperrectangle, p::Real=Inf)
diameter(H::Hyperrectangle, p::Real=Inf)
high(H::Hyperrectangle)
low(H::Hyperrectangle)
```

## VoidSets

```@docs
VoidSet
```

## Singletons

```@docs
Singleton
```

## Zonotopes

```@docs
Zonotope
Zonotope(center::Vector, G::AbstractMatrix)
dim(Z::Zonotope)
vertices_list(Z::Zonotope{Float64})
σ(d::AbstractVector, Z::Zonotope)
order(Z::Zonotope)
```
