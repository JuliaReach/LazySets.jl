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
dim(::Ball2{Float64})
σ(::AbstractVector{Float64}, ::Ball2{Float64})
∈(::AbstractVector{Float64}, ::Ball2{Float64})
an_element(::Ball2{Float64})
⊆(::Ball2{Float64}, ::Singleton{Float64})
⊆(::Ball2{Float64}, ::AbstractHyperrectangle{Float64})
is_intersection_empty(::Ball2{Float64}, ::Ball2{Float64}, ::Bool)
center(::Ball2{Float64})
```

### Infinity norm ball

```@docs
BallInf
dim(::BallInf{Float64})
σ(::AbstractVector{Float64}, ::BallInf{Float64})
∈(::AbstractVector{Float64}, ::BallInf{Float64})
an_element(::BallInf{Float64})
⊆(::BallInf{Float64}, ::AbstractHyperrectangle{Float64})
⊆(::BallInf{Float64}, ::LazySet{Float64})
is_intersection_empty(::BallInf{Float64}, ::AbstractHyperrectangle{Float64}, ::Bool)
norm(::BallInf{Float64}, ::Real)
radius(::BallInf{Float64}, ::Real)
diameter(::BallInf{Float64}, ::Real)
vertices_list(::BallInf{Float64})
singleton_list(::BallInf{Float64})
center(::BallInf{Float64})
radius_hyperrectangle(::BallInf{Float64})
radius_hyperrectangle(::BallInf{Float64}, ::Int)
```

### Manhattan norm ball

```@docs
Ball1
dim(::Ball1{Float64})
σ(::AbstractVector{Float64}, ::Ball1{Float64})
∈(::AbstractVector{Float64}, ::Ball1{Float64})
an_element(::Ball1{Float64})
⊆(::Ball1{Float64}, ::LazySet{Float64})
vertices_list(::Ball1{Float64})
singleton_list(::Ball1{Float64})
center(::Ball1{Float64})
```

### p-norm ball

```@docs
Ballp
dim(::Ballp{Float64})
σ(::AbstractVector{Float64}, ::Ballp{Float64})
∈(::AbstractVector{Float64}, ::Ballp{Float64})
an_element(::Ballp{Float64})
⊆(::Ballp{Float64}, ::Singleton{Float64})
⊆(::Ballp{Float64}, ::AbstractHyperrectangle{Float64})
center(::Ballp{Float64})
```

## Ellipsoid

```@docs
Ellipsoid
σ(::AbstractVector{Float64}, ::Ellipsoid{Float64})
center(::Ellipsoid{Float64})
∈(::AbstractVector{Float64}, ::Ellipsoid{Float64})
```

## Polygons

### Constraint representation

```@docs
HPolygon
dim(::HPolygon{Float64})
σ(::AbstractVector{Float64}, ::HPolygon{Float64})
∈(::AbstractVector{Float64}, ::HPolygon{Float64})
an_element(::HPolygon{Float64})
⊆(::HPolygon{Float64}, ::LazySet{Float64})
vertices_list(::HPolygon{Float64})
singleton_list(::HPolygon{Float64})
tohrep(::HPolygon{Float64})
tovrep(::HPolygon{Float64})
addconstraint!(::HPolygon{Float64}, ::LinearConstraint{Float64})
```

### Optimized constraint representation

```@docs
HPolygonOpt
dim(::HPolygonOpt{Float64})
σ(::AbstractVector{Float64}, ::HPolygonOpt{Float64})
∈(::AbstractVector{Float64}, ::HPolygonOpt{Float64})
an_element(::HPolygonOpt{Float64})
⊆(::HPolygonOpt{Float64}, ::LazySet{Float64})
vertices_list(::HPolygonOpt{Float64})
singleton_list(::HPolygonOpt{Float64})
tohrep(::HPolygonOpt{Float64})
tovrep(::HPolygonOpt{Float64})
addconstraint!(::HPolygonOpt{Float64}, ::LinearConstraint{Float64})
```

### Vertex representation

```@docs
VPolygon
dim(::VPolygon{Float64})
σ(::AbstractVector{Float64}, ::VPolygon{Float64})
∈(::AbstractVector{Float64}, ::VPolygon{Float64})
an_element(::VPolygon{Float64})
⊆(::VPolygon{Float64}, ::LazySet{Float64})
vertices_list(::VPolygon{Float64})
singleton_list(::VPolygon{Float64})
tohrep(::VPolygon{Float64})
tovrep(::VPolygon{Float64})
```

## Polytopes

```@docs
HPolytope
dim(P::HPolytope)
addconstraint!(P::HPolytope{Float64}, constraint::LinearConstraint{Float64})
constraints_list(P::HPolytope{Float64})
σ(d::AbstractVector{Float64}, P::HPolytope)
∈(::AbstractVector{Float64}, ::HPolytope{Float64})
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
dim(::Hyperrectangle{Float64})
σ(::AbstractVector{Float64}, ::Hyperrectangle{Float64})
∈(::AbstractVector{Float64}, ::Hyperrectangle{Float64})
an_element(::Hyperrectangle{Float64})
⊆(::Hyperrectangle{Float64}, ::AbstractHyperrectangle{Float64})
⊆(::Hyperrectangle{Float64}, ::LazySet{Float64})
is_intersection_empty(::Hyperrectangle{Float64}, ::AbstractHyperrectangle{Float64}, ::Bool)
norm(::Hyperrectangle{Float64}, ::Real)
radius(::Hyperrectangle{Float64}, ::Real)
diameter(::Hyperrectangle{Float64}, ::Real)
vertices_list(::Hyperrectangle{Float64})
singleton_list(::Hyperrectangle{Float64})
center(::Hyperrectangle{Float64})
radius_hyperrectangle(::Hyperrectangle{Float64})
radius_hyperrectangle(::Hyperrectangle{Float64}, ::Int)
high(::Hyperrectangle{Float64})
low(::Hyperrectangle{Float64})
```

## EmptySet

```@docs
EmptySet
dim(::EmptySet{Float64})
σ(::AbstractVector{Float64}, ::EmptySet{Float64})
∈(::AbstractVector{Float64}, ::EmptySet{Float64})
an_element(::EmptySet{Float64})
```

## Singletons

```@docs
Singleton
dim(::Singleton{Float64})
σ(::AbstractVector{Float64}, ::Singleton{Float64})
∈(::AbstractVector{Float64}, ::Singleton{Float64})
⊆(::Singleton{Float64}, ::AbstractSingleton{Float64})
⊆(::Singleton{Float64}, ::LazySet{Float64})
is_intersection_empty(::Singleton{Float64}, ::LazySet{Float64}, ::Bool)
is_intersection_empty(::LazySet{Float64}, ::Singleton{Float64}, ::Bool)
is_intersection_empty(::Singleton{Float64}, ::Singleton{Float64}, ::Bool)
norm(::Singleton{Float64}, ::Real)
diameter(::Singleton{Float64}, ::Real)
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
dim(::ZeroSet{Float64})
σ(::AbstractVector{Float64}, ::ZeroSet{Float64})
∈(::AbstractVector{Float64}, ::ZeroSet{Float64})
⊆(::ZeroSet{Float64}, ::AbstractSingleton{Float64})
⊆(::ZeroSet{Float64}, ::LazySet{Float64})
is_intersection_empty(::ZeroSet{Float64}, ::LazySet{Float64}, ::Bool)
is_intersection_empty(::LazySet{Float64}, ::ZeroSet{Float64}, ::Bool)
is_intersection_empty(::ZeroSet{Float64}, ::ZeroSet{Float64}, ::Bool)
norm(::ZeroSet{Float64}, ::Real)
diameter(::ZeroSet{Float64}, ::Real)
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
dim(::Zonotope{Float64})
σ(::AbstractVector{Float64}, ::Zonotope{Float64})
∈(::AbstractVector{Float64}, ::Zonotope{Float64})
an_element(::Zonotope{Float64})
⊆(::Zonotope{Float64}, ::LazySet{Float64})
center(::Zonotope{Float64})
vertices_list(::Zonotope{Float64})
singleton_list(::Zonotope{Float64})
order(::Zonotope{Float64})
minkowski_sum(Z1::Zonotope, Z2::Zonotope)
linear_map(M::AbstractMatrix, Z::Zonotope)
scale(α::Real, Z::Zonotope)
```

## Hyperplane

```@docs
Hyperplane
dim(::Hyperplane{Float64})
σ(::AbstractVector{Float64}, ::Hyperplane{Float64})
an_element(::Hyperplane{Float64})
∈(::AbstractVector{Float64}, ::Hyperplane{Float64})
```

## Half-Space

```@docs
HalfSpace
dim(::HalfSpace{Float64})
σ(::AbstractVector{Float64}, ::HalfSpace{Float64})
an_element(::HalfSpace{Float64})
∈(::AbstractVector{Float64}, ::HalfSpace{Float64})
```
