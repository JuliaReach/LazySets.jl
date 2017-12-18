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

## Support function and support vector

```@docs
LazySets
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
an_element(::Ball2{Float64})
center(::Ball2{Float64})
```

### Infinity norm ball

```@docs
BallInf
dim(::BallInf)
σ(::AbstractVector{Float64}, ::BallInf{Float64})
∈(::AbstractVector{Float64}, ::BallInf{Float64})
an_element(::BallInf{Float64})
⊆(::BallInf, ::AbstractHyperrectangle)
norm(::BallInf, ::Real=Inf)
radius(::BallInf, ::Real=Inf)
diameter(::BallInf, ::Real=Inf)
vertices_list(::BallInf{Float64})
singleton_list(::BallInf{Float64})
center(::BallInf{Float64})
radius_hyperrectangle(::BallInf{Float64})
radius_hyperrectangle(::BallInf{Float64}, ::Int)
```

### Manhattan norm ball

```@docs
Ball1
dim(::Ball1)
σ(::AbstractVector{Float64}, ::Ball1{Float64})
∈(::AbstractVector{Float64}, ::Ball1{Float64})
an_element(::Ball1{Float64})
vertices_list(::Ball1{Float64})
singleton_list(::Ball1{Float64})
center(::Ball1{Float64})
```

### p-norm ball

```@docs
Ballp
dim(::Ballp)
σ(::AbstractVector{Float64}, ::Ballp)
∈(::AbstractVector{Float64}, ::Ballp{Float64})
an_element(::Ballp{Float64})
center(::Ballp{Float64})
```

## Polygons

### Constraint representation

```@docs
HPolygon
dim(::HPolygon)
σ(::AbstractVector{Float64}, ::HPolygon{Float64})
∈(::AbstractVector{Float64}, ::HPolygon{Float64})
an_element(::HPolygon{Float64})
vertices_list(::HPolygon{Float64})
singleton_list(::HPolygon{Float64})
tohrep(::HPolygon{Float64})
tovrep(::HPolygon{Float64})
addconstraint!(::HPolygon{Float64}, ::LinearConstraint{Float64})
```

### Optimized constraint representation

```@docs
HPolygonOpt
dim(::HPolygonOpt)
σ(::AbstractVector{Float64}, ::HPolygonOpt{Float64})
∈(::AbstractVector{Float64}, ::HPolygonOpt{Float64})
an_element(::HPolygonOpt{Float64})
vertices_list(::HPolygonOpt{Float64})
singleton_list(::HPolygonOpt{Float64})
tohrep(::HPolygonOpt{Float64})
tovrep(::HPolygonOpt{Float64})
addconstraint!(::HPolygonOpt{Float64}, ::LinearConstraint{Float64})
```

### Vertex representation

```@docs
VPolygon
dim(::VPolygon)
σ(::AbstractVector{Float64}, ::VPolygon{Float64})
∈(::AbstractVector{Float64}, ::VPolygon{Float64})
an_element(::VPolygon{Float64})
vertices_list(::VPolygon{Float64})
singleton_list(::VPolygon{Float64})
tohrep(::VPolygon{Float64})
tovrep(::VPolygon{Float64})
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
σ(::AbstractVector{Float64}, ::Hyperrectangle{Float64})
∈(::AbstractVector{Float64}, ::Hyperrectangle{Float64})
an_element(::Hyperrectangle{Float64})
⊆(::Hyperrectangle, ::AbstractHyperrectangle)
norm(::Hyperrectangle, ::Real=Inf)
radius(::Hyperrectangle, ::Real=Inf)
diameter(::Hyperrectangle, ::Real=Inf)
vertices_list(::Hyperrectangle{Float64})
singleton_list(::Hyperrectangle{Float64})
center(::Hyperrectangle{Float64})
radius_hyperrectangle(::Hyperrectangle{Float64})
radius_hyperrectangle(::Hyperrectangle{Float64}, ::Int)
high(::Hyperrectangle)
low(::Hyperrectangle)
```

## EmptySet

```@docs
EmptySet
dim(::EmptySet)
σ(::AbstractVector{Float64}, ::EmptySet)
∈(::AbstractVector{Float64}, ::EmptySet)
an_element(::EmptySet)
```

## Singletons

```@docs
Singleton
dim(::Singleton)
σ(::AbstractVector{Float64}, ::Singleton{Float64})
∈(::AbstractVector{Float64}, ::Singleton{Float64})
⊆(::Singleton, ::LazySet)
norm(::Singleton, ::Real=Inf)
diameter(::Singleton, ::Real=Inf)
vertices_list(::Singleton{Float64})
singleton_list(::Singleton{Float64})
center(::Singleton{Float64})
radius_hyperrectangle(::Singleton{Float64})
radius_hyperrectangle(::Singleton{Float64}, ::Int)
an_element(::Singleton{Float64})
element(::Singleton{Float64})
element(::Singleton{Float64}, ::Int)
```

### ZeroSet

```@docs
ZeroSet
dim(::ZeroSet)
σ(::AbstractVector{Float64}, ::ZeroSet)
∈(::AbstractVector{Float64}, ::ZeroSet{Float64})
⊆(::ZeroSet, ::LazySet)
norm(::ZeroSet, ::Real=Inf)
diameter(::ZeroSet, ::Real=Inf)
vertices_list(::ZeroSet{Float64})
singleton_list(::ZeroSet{Float64})
center(::ZeroSet{Float64})
radius_hyperrectangle(::ZeroSet{Float64})
radius_hyperrectangle(::ZeroSet{Float64}, ::Int)
an_element(::ZeroSet{Float64})
element(::ZeroSet{Float64})
element(::ZeroSet{Float64}, ::Int)
```

## Zonotopes

```@docs
Zonotope
dim(::Zonotope)
σ(::AbstractVector{Float64}, Z::Zonotope)
∈(::AbstractVector{Float64}, ::Zonotope{Float64})
an_element(::Zonotope{Float64})
center(::Zonotope{Float64})
vertices_list(::Zonotope{Float64})
singleton_list(::Zonotope{Float64})
order(::Zonotope)
```
