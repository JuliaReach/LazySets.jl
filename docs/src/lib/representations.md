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
dim(::Ball2)
σ(::AbstractVector{Float64}, ::Ball2)
```

### Infinity norm ball

```@docs
BallInf
dim(::BallInf)
σ(::AbstractVector{Float64}, ::BallInf)
vertices_list(::BallInf)
norm(::BallInf, ::Real=Inf)
radius(::BallInf, ::Real=Inf)
diameter(::BallInf, ::Real=Inf)
```

## Polygons

### Constraint representation

```@docs
HPolygon
addconstraint!(::HPolygon{Float64}, ::LinearConstraint{Float64})
dim(::HPolygon)
σ(::AbstractVector{Float64}, ::HPolygon)
is_contained(::AbstractVector{Float64}, ::HPolygon{Float64})
tovrep(::HPolygon)
vertices_list(::HPolygon)
```

### Optimized constraint representation

```@docs
HPolygonOpt
addconstraint!(::HPolygonOpt{Float64}, ::LinearConstraint{Float64})
dim(::HPolygonOpt)
σ(::AbstractVector{Float64}, ::HPolygonOpt{Float64})
is_contained(::AbstractVector{Float64}, ::HPolygonOpt)
tovrep(::HPolygonOpt)
vertices_list(::HPolygonOpt)
```

### Vertex representation

```@docs
VPolygon
dim(::VPolygon)
σ(::AbstractVector{Float64}, ::VPolygon)
vertices_list(::VPolygon)
singleton_list(::VPolygon)
```

## Lines and linear constraints

```@docs
LinearConstraint
Line
intersection(::Line{Float64}, L2::Line{Float64})
```

## Hyperrectangles

```@docs
Hyperrectangle
Hyperrectangle(;kwargs...)
dim(::Hyperrectangle)
σ(::AbstractVector{Float64}, ::Hyperrectangle)
vertices_list(::Hyperrectangle)
norm(::Hyperrectangle, ::Real=Inf)
radius(::Hyperrectangle, ::Real=Inf)
diameter(::Hyperrectangle, ::Real=Inf)
high(::Hyperrectangle)
low(::Hyperrectangle)
```

## VoidSets

```@docs
VoidSet
```

## Singletons

```@docs
Singleton
dim(::Singleton)
σ(::AbstractVector{Float64}, ::Singleton)
```

## Zonotopes

```@docs
Zonotope
dim(::Zonotope)
σ(d::AbstractVector{Float64}, Z::Zonotope)
vertices_list(::Zonotope{Float64})
order(::Zonotope)
```
