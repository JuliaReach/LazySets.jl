# Common Set Representations

This section of the manual describes the basic set representation types.

```@contents
Pages = ["representations.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

## Balls

### Euclidean norm ball

```@docs
Ball2
σ(::AbstractVector{N}, ::Ball2{N}) where {N<:AbstractFloat}
∈(::AbstractVector{N}, ::Ball2{N}) where {N<:AbstractFloat}
center(::Ball2{N}) where {N<:AbstractFloat}
rand(::Type{Ball2})
sample(::Ball2{N}, ::Int) where {N<:AbstractFloat}
translate(::Ball2{N}, ::AbstractVector{N}) where {N<:AbstractFloat}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`isbounded`](@ref isbounded(::AbstractCentrallySymmetric))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetric{N}) where {N<:Real})

### Infinity norm ball

```@docs
BallInf
center(::BallInf{N}) where {N<:Real}
radius(::BallInf, ::Real=Inf)
radius_hyperrectangle(::BallInf{N}) where {N<:Real}
radius_hyperrectangle(::BallInf{N}, ::Int) where {N<:Real}
isflat(::BallInf)
rand(::Type{BallInf})
σ(::AbstractVector{N}, ::BallInf{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::BallInf{N}) where {N<:Real}
translate(::BallInf{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope{N}) where {N<:Real})

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`order`](@ref order(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle{N}) where {N<:Real})
* [`high`](@ref high(::AbstractHyperrectangle{N}) where {N<:Real})
* [`low`](@ref low(::AbstractHyperrectangle{N}) where {N<:Real})
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N<:Real})

### Manhattan norm ball

```@docs
Ball1
σ(::AbstractVector{N}, ::Ball1{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Ball1{N}) where {N<:Real}
vertices_list(::Ball1{N}) where {N<:Real}
center(::Ball1{N}) where {N<:Real}
rand(::Type{Ball1})
constraints_list(::Ball1{N}) where {N<:Real}
translate(::Ball1{N}, ::AbstractVector{N}) where {N<:Real}
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope{N}) where {N<:Real})

### p-norm ball

```@docs
Ballp
σ(::AbstractVector{N}, ::Ballp{N}) where {N<:AbstractFloat}
∈(::AbstractVector{N}, ::Ballp{N}) where {N<:AbstractFloat}
center(::Ballp{N}) where {N<:AbstractFloat}
rand(::Type{Ballp})
translate(::Ballp{N}, ::AbstractVector{N}) where {N<:AbstractFloat}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`isbounded`](@ref isbounded(::AbstractCentrallySymmetric))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetric{N}) where {N<:Real})

## Ellipsoid

```@docs
Ellipsoid
ρ(::AbstractVector{N}, ::Ellipsoid{N}) where {N<:AbstractFloat}
σ(::AbstractVector{N}, ::Ellipsoid{N}) where {N<:AbstractFloat}
∈(::AbstractVector{N}, ::Ellipsoid{N}) where {N<:AbstractFloat}
rand(::Type{Ellipsoid})
center(::Ellipsoid{N}) where {N<:AbstractFloat}
translate(::Ellipsoid{N}, ::AbstractVector{N}) where {N<:AbstractFloat}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractCentrallySymmetric`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetric))
* [`isbounded`](@ref isbounded(::AbstractCentrallySymmetric))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetric))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetric{N}) where {N<:Real})

## Empty set

```@docs
EmptySet
∅
dim(::EmptySet)
σ(::AbstractVector{N}, ::EmptySet{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::EmptySet{N}) where {N<:Real}
∈(::AbstractVector{N}, ::EmptySet{N}) where {N<:Real}
an_element(::EmptySet)
rand(::Type{EmptySet})
isbounded(::EmptySet)
isempty(::EmptySet)
norm(::EmptySet, ::Real=Inf)
radius(::EmptySet, ::Real=Inf)
diameter(::EmptySet, ::Real=Inf)
linear_map(::AbstractMatrix{N}, ::EmptySet{N}) where {N}
translate(::EmptySet{N}, ::AbstractVector{N}) where {N<:Real}
plot_recipe(::EmptySet{N}, ::N=zero(N)) where {N<:Real}
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::EmptySet{N}, ::N=zero(N)) where {N<:Real}
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
ρ(::AbstractVector{N}, ::HalfSpace{N}) where {N<:Real}
σ(::AbstractVector{N}, ::HalfSpace{N}) where {N<:Real}
∈(::AbstractVector{N}, ::HalfSpace{N}) where {N<:Real}
an_element(::HalfSpace{N}) where {N<:Real}
rand(::Type{HalfSpace})
isbounded(::HalfSpace)
isuniversal(::HalfSpace{N}, ::Bool=false) where {N<:Real}
isempty(::HalfSpace)
constraints_list(::HalfSpace{N}) where {N<:Real}
constraints_list(::AbstractMatrix{N}, ::AbstractVector{N}) where {N<:Real}
constrained_dimensions(::HalfSpace{N}) where {N<:Real}
translate(::HalfSpace{N}, ::AbstractVector{N}) where {N<:Real}
halfspace_left(::AbstractVector{N}, ::AbstractVector{N}) where {N<:Real}
halfspace_right(::AbstractVector{N}, ::AbstractVector{N}) where {N<:Real}
tosimplehrep(::AbstractVector{LC}) where {N<:Real, LC<:LinearConstraint{N}}
remove_redundant_constraints(::AbstractVector{LC}) where {N<:Real, LC<:LinearConstraint{N}}
remove_redundant_constraints!(::AbstractVector{LC}) where {N<:Real, LC<:LinearConstraint{N}}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

## Hyperplane

```@docs
Hyperplane
dim(::Hyperplane)
ρ(::AbstractVector{N}, ::Hyperplane{N}) where {N<:Real}
σ(::AbstractVector{N}, ::Hyperplane{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Hyperplane{N}) where {N<:Real}
an_element(::Hyperplane{N}) where {N<:Real}
rand(::Type{Hyperplane})
isbounded(::Hyperplane)
isuniversal(::Hyperplane{N}, ::Bool=false) where {N<:Real}
isempty(::Hyperplane)
constrained_dimensions(::Hyperplane{N}) where {N<:Real}
constraints_list(::Hyperplane{N}) where {N<:Real}
translate(::Hyperplane{N}, ::AbstractVector{N}) where {N<:Real}
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
translate(::Hyperrectangle{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope{N}) where {N<:Real})

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`order`](@ref order(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`σ`](@ref σ(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`ρ`](@ref ρ(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractHyperrectangle{N}) where {N<:Real})
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`vertices_list`](@ref vertices_list(::AbstractHyperrectangle{N}) where {N<:Real})
* [`high`](@ref high(::AbstractHyperrectangle{N}) where {N<:Real})
* [`low`](@ref low(::AbstractHyperrectangle{N}) where {N<:Real})
* [`isflat`](@ref isflat(::Hyperrectangle))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N<:Real})

## Interval

```@docs
Interval
dim(::Interval)
σ(::AbstractVector{N}, ::Interval{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::Interval{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Interval{N}) where {N<:Real}
∈(::N, ::Interval{N}) where {N<:Real}
an_element(::Interval{N}) where {N<:Real}
vertices_list(::Interval{N}) where {N<:Real}
translate(::Interval{N}, ::AbstractVector{N}) where {N<:Real}
center(::Interval{N}) where {N<:Real}
min(::Interval{N}) where {N<:Real}
max(::Interval{N}) where {N<:Real}
low(::Interval{N}) where {N<:Real}
high(::Interval{N}) where {N<:Real}
radius_hyperrectangle(::Interval{N}) where {N<:Real}
radius_hyperrectangle(::Interval{N}, ::Int) where {N<:Real}
+(::Interval{N}, ::Interval{N}) where {N<:Real}
-(::Interval{N}, ::Interval{N}) where {N<:Real}
*(::Interval{N}, ::Interval{N}) where {N<:Real}
rand(::Type{Interval})
isflat(::Interval)
plot_recipe(::Interval{N}, ::N=zero(N)) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`order`](@ref order(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`generators`](@ref generators(::AbstractHyperrectangle))
* [`genmat`](@ref genmat(::AbstractHyperrectangle))
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N<:Real})

## Line

```@docs
Line
dim(::Line)
σ(::AbstractVector{N}, ::Line{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Line{N}) where {N<:Real}
an_element(::Line{N}) where {N<:Real}
rand(::Type{Line})
isbounded(::Line)
isuniversal(::Line{N}, ::Bool=false) where {N<:Real}
isempty(::Line)
constrained_dimensions(::Line{N}) where {N<:Real}
constraints_list(::Line{N}) where {N<:Real}
translate(::Line{N}, ::AbstractVector{N}) where {N<:Real}
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
translate(::LineSegment{N}, ::AbstractVector{N}) where {N<:Real}
generators(::LineSegment{N}) where {N<:Real}
genmat(::LineSegment)
plot_recipe(::LineSegment{N}, ::N=zero(N)) where {N<:Real}
RecipesBase.apply_recipe(::Dict{Symbol,Any}, ::Union{LineSegment{N}, Interval{N}}, ::N=zero(N)) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`order`](@ref order(::AbstractZonotope))

## Polygons

### Constraint representation

```@docs
HPolygon
σ(::AbstractVector{N}, ::HPolygon{N}) where {N<:Real}
translate(::HPolygon{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractPolyhedron{N}) where {N<:Real})

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))

Inherited from [`AbstractHPolygon`](@ref):
* [`an_element`](@ref an_element(::AbstractHPolygon{N}) where {N<:Real})
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractHPolygon{N}) where {N<:Real})
* [`vertices_list`](@ref vertices_list(::AbstractHPolygon{N}, ::Bool=false, ::Bool=true) where {N<:Real})
* [`tohrep`](@ref tohrep(::HPOLYGON) where {HPOLYGON<:AbstractHPolygon})
* [`tovrep`](@ref tovrep(::AbstractHPolygon{N}) where {N<:Real})
* [`isbounded`](@ref isbounded(::AbstractHPolygon, ::Bool=true))
* [`addconstraint!`](@ref addconstraint!(::AbstractHPolygon{N}, ::LinearConstraint{N}) where {N<:Real})
* [`isredundant`](@ref isredundant(::LinearConstraint{N}, ::LinearConstraint{N}, ::LinearConstraint{N}) where {N<:Real})
* [`remove_redundant_constraints!`](@ref remove_redundant_constraints!(::AbstractHPolygon))
* [`constraints_list`](@ref constraints_list(::AbstractHPolygon{N}) where {N<:Real})

### Optimized constraint representation

```@docs
HPolygonOpt
σ(::AbstractVector{N}, ::HPolygonOpt{N}) where {N<:Real}
translate(::HPolygonOpt{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractPolygon`](@ref):
* [`dim`](@ref dim(::AbstractPolygon))

Inherited from [`AbstractHPolygon`](@ref):
* [`an_element`](@ref an_element(::AbstractHPolygon{N}) where {N<:Real})
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractHPolygon{N}) where {N<:Real})
* [`vertices_list`](@ref vertices_list(::AbstractHPolygon{N}, ::Bool=false, ::Bool=true) where {N<:Real})
* [`tohrep`](@ref tohrep(::HPOLYGON) where {HPOLYGON<:AbstractHPolygon})
* [`tovrep`](@ref tovrep(::AbstractHPolygon{N}) where {N<:Real})
* [`isbounded`](@ref isbounded(::AbstractHPolygon, ::Bool=true))
* [`addconstraint!`](@ref addconstraint!(::AbstractHPolygon{N}, ::LinearConstraint{N}) where {N<:Real})
* [`isredundant`](@ref isredundant(::LinearConstraint{N}, ::LinearConstraint{N}, ::LinearConstraint{N}) where {N<:Real})
* [`remove_redundant_constraints!`](@ref remove_redundant_constraints!(::AbstractHPolygon))
* [`constraints_list`](@ref constraints_list(::AbstractHPolygon{N}) where {N<:Real})

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
translate(::VPolygon{N}, ::AbstractVector{N}) where {N<:Real}
remove_redundant_vertices(::VPolygon{N}; ::String="monotone_chain") where {N<:Real}
remove_redundant_vertices!(::VPolygon{N}; ::String="monotone_chain") where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractPolyhedron{N}) where {N<:Real})

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
addconstraint!(::HPoly{N}, ::LinearConstraint{N}) where {N<:Real}
constraints_list(::HPoly{N}) where {N<:Real}
tohrep(::HPoly{N}) where {N<:Real}
tovrep(::HPoly{N}) where {N<:Real}
isempty(::HPoly{N}, ::Bool=false) where {N<:Real}
translate(::PT, ::AbstractVector{N}) where {N<:Real, PT<:HPoly{N}}
polyhedron(::HPoly{N}) where {N<:Real}
remove_redundant_constraints(::PT) where {N<:Real, PT<:HPoly{N}}
remove_redundant_constraints!(::HPoly{N}) where {N<:Real}
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolyhedron`](@ref):
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractPolyhedron{N}) where {N<:Real})
* [`constrained_dimensions`](@ref constrained_dimensions(::AbstractPolyhedron)
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractPolyhedron{N}) where {N<:Real})

#### Polytopes

The following methods are specific for `HPolytope`.

```@docs
rand(::Type{HPolytope})
vertices_list(::HPolytope{N}) where {N<:Real}
isbounded(::HPolytope, ::Bool=true)
```

Inherited from [`AbstractPolytope`](@ref):
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

#### Polyhedra

The following methods are specific for `HPolyhedron`.

```@docs
rand(::Type{HPolyhedron})
isbounded(::HPolyhedron)
isuniversal(::HPolyhedron{N}, ::Bool=false) where {N<:Real}
vertices_list(::HPolyhedron{N}) where {N<:Real}
singleton_list(::HPolyhedron{N}) where {N<:Real}
```

### Vertex representation

```@docs
VPolytope
dim(::VPolytope)
σ(::AbstractVector{N}, ::VPolytope{N}) where {N<:Real}
∈(::AbstractVector{N}, ::VPolytope{N}) where {N<:Real}
rand(::Type{VPolytope})
translate(::VPolytope{N}, ::AbstractVector{N}) where {N<:Real}
vertices_list(::VPolytope{N}) where {N<:Real}
remove_redundant_vertices(::VPolytope{N}) where {N<:Real}
constraints_list(::VPolytope{N}) where {N<:Real}
tohrep(::VPolytope{N}) where {N<:Real}
tovrep(::VPolytope)
polyhedron(::VPolytope{N}) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::VPolytope{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`isempty`](@ref isempty(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractPolyhedron{N}) where {N<:Real})

## Polynomial Zonotopes

```@docs
PolynomialZonotope
dim(::PolynomialZonotope)
σ(::AbstractVector{N}, ::PolynomialZonotope{N}) where {N}
ρ(::AbstractVector{N}, ::PolynomialZonotope{N}) where {N}
polynomial_order(pz::PolynomialZonotope)
order(::PolynomialZonotope)
linear_map(::Matrix, ::PolynomialZonotope)
scale(::Number, ::PolynomialZonotope)
```

## Singleton

```@docs
Singleton
rand(::Type{Singleton})
element(::Singleton{N}) where {N<:Real}
element(::Singleton{N}, ::Int) where {N<:Real}
translate(::Singleton{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`order`](@ref order(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`high`](@ref high(::AbstractHyperrectangle{N}) where {N<:Real})
* [`low`](@ref low(::AbstractHyperrectangle{N}) where {N<:Real})
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N<:Real})

Inherited from [`AbstractSingleton`](@ref):
* [`σ`](@ref σ(::AbstractVector{N}, ::AbstractSingleton{N}) where {N<:Real})
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractSingleton{N}) where {N<:Real})
* [`an_element`](@ref an_element(::AbstractSingleton{N}) where {N<:Real})
* [`center`](@ref center(::AbstractSingleton{N}) where {N<:Real})
* [`vertices_list`](@ref vertices_list(::AbstractSingleton{N}) where {N<:Real})
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}) where {N<:Real})
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}, ::Int) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractSingleton{N}) where {N<:Real})
* [`generators`](@ref generators(::AbstractSingleton{N}) where {N<:Real})
* [`genmat`](@ref genmat(::AbstractSingleton{N}) where {N<:Real})

## Universe

```@docs
Universe
dim(::Universe)
ρ(::AbstractVector{N}, ::Universe{N}) where {N<:Real}
σ(::AbstractVector{N}, ::Universe{N}) where {N<:Real}
∈(::AbstractVector{N}, ::Universe{N}) where {N<:Real}
rand(::Type{Universe})
an_element(::Universe{N}) where {N<:Real}
isempty(::Universe)
isbounded(::Universe)
isuniversal(::Universe{N}, ::Bool=false) where {N<:Real}
norm(::Universe, ::Real=Inf)
radius(::Universe, ::Real=Inf)
diameter(::Universe, ::Real=Inf)
constraints_list(::Universe{N}) where {N<:Real}
constrained_dimensions(::Universe)
translate(::Universe{N}, ::AbstractVector{N}) where {N<:Real}
```

## Zero set

```@docs
ZeroSet
dim(::ZeroSet)
σ(::AbstractVector{N}, ::ZeroSet{N}) where {N<:Real}
ρ(::AbstractVector{N}, ::ZeroSet{N}) where {N<:Real}
∈(::AbstractVector{N}, ::ZeroSet{N}) where {N<:Real}
rand(::Type{ZeroSet})
element(::ZeroSet{N}) where {N<:Real}
element(::ZeroSet{N}, ::Int) where {N<:Real}
linear_map(::AbstractMatrix{N}, ::ZeroSet{N}) where {N<:Real}
translate(::ZeroSet{N}, ::AbstractVector{N}) where {N<:Real}
```
Inherited from [`LazySet`](@ref):
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))

Inherited from [`AbstractZonotope`](@ref):
* [`ngens`](@ref ngens(::AbstractZonotope))
* [`order`](@ref order(::AbstractZonotope))

Inherited from [`AbstractHyperrectangle`](@ref):
* [`norm`](@ref norm(::AbstractHyperrectangle, ::Real))
* [`radius`](@ref radius(::AbstractHyperrectangle, ::Real))
* [`high`](@ref high(::AbstractHyperrectangle{N}) where {N<:Real})
* [`low`](@ref low(::AbstractHyperrectangle{N}) where {N<:Real})
* [`constraints_list`](@ref constraints_list(::AbstractHyperrectangle{N}) where {N<:Real})

Inherited from [`AbstractSingleton`](@ref):
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}) where {N<:Real})
* [`radius_hyperrectangle`](@ref radius_hyperrectangle(::AbstractSingleton{N}, ::Int) where {N<:Real})
* [`vertices_list`](@ref vertices_list(::AbstractSingleton{N}) where {N<:Real})
* [`center`](@ref center(::AbstractSingleton{N}) where {N<:Real})
* [`an_element`](@ref an_element(::AbstractSingleton{N}) where {N<:Real})
* [`generators`](@ref generators(::AbstractSingleton{N}) where {N<:Real})
* [`genmat`](@ref genmat(::AbstractSingleton{N}) where {N<:Real})

## Zonotope

```@docs
Zonotope
center(::Zonotope{N}) where {N<:Real}
rand(::Type{Zonotope})
generators(Z::Zonotope)
genmat(Z::Zonotope)
scale(::Real, ::Zonotope)
ngens(::Zonotope)
reduce_order(::Zonotope, r)
split(::Zonotope, ::Int)
```

Inherited from [`LazySet`](@ref):
* [`norm`](@ref norm(::LazySet, ::Real))
* [`radius`](@ref radius(::LazySet, ::Real))
* [`diameter`](@ref diameter(::LazySet, ::Real))

Inherited from [`AbstractPolytope`](@ref):
* [`isbounded`](@ref isbounded(::AbstractPolytope))
* [`singleton_list`](@ref singleton_list(::AbstractPolytope{N}) where {N<:Real})

Inherited from [`AbstractCentrallySymmetricPolytope`](@ref):
* [`dim`](@ref dim(::AbstractCentrallySymmetricPolytope))
* [`isempty`](@ref isempty(::AbstractCentrallySymmetricPolytope))
* [`an_element`](@ref an_element(::AbstractCentrallySymmetricPolytope{N}) where {N<:Real})

Inherited from [`AbstractZonotope`](@ref):
* [`ρ`](@ref ρ(::AbstractVector{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`σ`](@ref σ(::AbstractVector{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`∈`](@ref ∈(::AbstractVector{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`linear_map`](@ref linear_map(::AbstractMatrix{N}, ::AbstractZonotope{N}) where {N<:Real})
* [`translate`](@ref translate(::AbstractZonotope{N}, ::AbstractVector{N}) where {N<:Real})
* [`constraints_list`](@ref constraints_list(::AbstractZonotope{N}) where {N<:Real})
* [`constraints_list`](@ref constraints_list(::AbstractZonotope{N}; ::Bool=true) where {N<:AbstractFloat})
* [`vertices_list`](@ref vertices_list(::AbstractZonotope{N}) where {N<:Real})
* [`order`](@ref order(::AbstractZonotope))
