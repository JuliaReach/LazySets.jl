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

Unit balls are defined by int center (vector) and radius (scalar), such as
infinity-norm balls,

$B_\infty(c, r) = \{ x ∈ \mathbb{R}^n : \Vert x - c\Vert_\infty \leq r \}.$

and Euclidean (2-norm) balls,

$B_2(c, r) = \{ x ∈ \mathbb{R}^n : \Vert x - c\Vert_2 \leq r \}.$

```@docs
BallInf
dim(B::BallInf)
σ(d::AbstractVector{Float64}, B::BallInf)
Ball2
dim(B::Ball2)
σ(d::AbstractVector{Float64}, B::Ball2)
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
radius(H::Hyperrectangle)
dim(H::Hyperrectangle)
σ(d::AbstractVector{Float64}, H::Hyperrectangle)
vertices_list(H::Hyperrectangle)
diameter(H::Hyperrectangle)
```

## VoidSets

```@docs
VoidSet
```

## Singletons

```@docs
Singleton
```
