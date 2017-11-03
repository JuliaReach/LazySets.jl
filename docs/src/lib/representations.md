# Common Set Representations

This section of the manual describes the basic set representation types.

```@contents
Pages = ["representations.md"]
```

```@meta
CurrentModule = LazySets
```

```@docs
LazySets
LazySets.LazySet
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
tovrep(s::HPolygon)
addconstraint!(p::HPolygon, c::LinearConstraint)
σ(d::AbstractVector{Float64}, p::HPolygon)
```

### Vertex representation

```@docs
VPolygon
```

### Optimized constraint representation

```@docs
HPolygonOpt
σ(d::AbstractVector{Float64}, p::HPolygonOpt)
tovrep(po::HPolygonOpt)
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
```

## VoidSets

```@docs
VoidSet
```

## Singletons

```@docs
Singleton
```
