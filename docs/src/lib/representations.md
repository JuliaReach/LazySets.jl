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

## Abstract support function and support vector

```@docs
LazySets
LazySets.LazySet
ρ
support_function
support_vector
```

## Balls

### Euclidean norm ball

```@docs
Ball2
dim(::Ball2)
σ(::AbstractVector{Float64}, ::Ball2)
∈(::AbstractVector{Float64}, ::Ball2{Float64})
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
∈(::AbstractVector{Float64}, ::BallInf{Float64})
```

### Manhattan norm ball

```@docs
Ball1
dim(::Ball1)
σ(::AbstractVector{Float64}, ::Ball1)
∈(::AbstractVector{Float64}, ::Ball1{Float64})
```

### p-norm ball

```@docs
Ballp
dim(::Ballp)
σ(::AbstractVector{Float64}, ::Ballp)
∈(::AbstractVector{Float64}, ::Ballp{Float64})
```

## Polygons

### Constraint representation

```@docs
HPolygon
addconstraint!(::HPolygon{Float64}, ::LinearConstraint{Float64})
dim(::HPolygon)
σ(::AbstractVector{Float64}, ::HPolygon)
∈(::AbstractVector{Float64}, ::HPolygon{Float64})
tovrep(::HPolygon)
vertices_list(::HPolygon)
```

### Optimized constraint representation

```@docs
HPolygonOpt
addconstraint!(::HPolygonOpt{Float64}, ::LinearConstraint{Float64})
dim(::HPolygonOpt)
σ(::AbstractVector{Float64}, ::HPolygonOpt{Float64})
∈(::AbstractVector{Float64}, ::HPolygonOpt{Float64})
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
∈(::AbstractVector{Float64}, ::VPolygon{Float64})
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
∈(::AbstractVector{Float64}, ::Hyperrectangle{Float64})
high(::Hyperrectangle)
low(::Hyperrectangle)
```

## EmptySet

```@docs
EmptySet
dim(::EmptySet)
σ(::AbstractVector{Float64}, ::EmptySet)
∈(::AbstractVector{Float64}, ::EmptySet)
```

## ZeroSet

```@docs
ZeroSet
dim(::ZeroSet)
σ(::AbstractVector{Float64}, ::ZeroSet)
∈(::AbstractVector{Float64}, ::ZeroSet)
```

## Singletons

```@docs
Singleton
dim(::Singleton)
σ(::AbstractVector{Float64}, ::Singleton)
∈(::AbstractVector{Float64}, ::Singleton{Float64})
⊆(::Singleton, ::LazySet)
```

## Zonotopes

```@docs
Zonotope
dim(::Zonotope)
σ(d::AbstractVector{Float64}, Z::Zonotope)
vertices_list(::Zonotope{Float64})
order(::Zonotope)
∈(::AbstractVector{Float64}, ::Zonotope{Float64})
```
