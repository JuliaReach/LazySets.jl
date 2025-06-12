# Conversion between set representations

This section of the manual lists the conversion functions between set
representations.

```@contents
Pages = ["conversion.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

```@docs
convert(::Type{Interval}, ::Rectification{N, IN}) where {N, IN<:Interval}
convert(::Type{Interval}, ::MinkowskiSum{N, IT, IT}) where {N, IT<:Interval}
convert(::Type{Hyperrectangle}, ::CartesianProduct{N, HN1, HN2}) where {N, HN1<:AbstractHyperrectangle, HN2<:AbstractHyperrectangle}
convert(::Type{Hyperrectangle}, ::CartesianProductArray{N, HN}) where {N, HN<:AbstractHyperrectangle}
convert(::Type{Hyperrectangle}, ::CartesianProductArray{N, IN}) where {N, IN<:Interval}
convert(::Type{Hyperrectangle}, ::Rectification{N, AH}) where {N, AH<:AbstractHyperrectangle}
convert(::Type{HPolygon}, ::LazySet)
convert(::Type{HPolygon}, ::VPolygon)
convert(::Type{HPolygon}, ::LineSegment{N}) where {N}
convert(::Type{HPolygonOpt}, ::LazySet)
convert(::Type{HPolygonOpt}, ::VPolygon)
convert(::Type{HPolygonOpt}, ::LineSegment{N}) where {N}
convert(::Type{VPolygon}, ::AbstractHPolygon)
convert(::Type{Zonotope}, ::LinearMap{N, ZN}) where {N, ZN<:AbstractZonotope}
convert(::Type{Zonotope}, ::LinearMap{N, CartesianProduct{N, HN1, HN2}}) where {N, HN1<:AbstractHyperrectangle, HN2<:AbstractHyperrectangle}
convert(::Type{Zonotope}, ::LinearMap{N, CartesianProductArray{N, HN}}) where {N, HN<:AbstractHyperrectangle}
convert(::Type{Zonotope}, ::CartesianProduct{N, ZN1, ZN2}) where {N, ZN1<:AbstractZonotope, ZN2<:AbstractZonotope}
convert(::Type{Zonotope}, ::CartesianProduct{N, HN1, HN2}) where {N, HN1<:AbstractHyperrectangle, HN2<:AbstractHyperrectangle}
convert(::Type{Zonotope}, ::CartesianProductArray{N, AZ}) where {N, AZ<:AbstractZonotope}
convert(::Type{Zonotope}, ::CartesianProductArray{N, HN}) where {N, HN<:AbstractHyperrectangle}
convert(::Type{CartesianProduct{N, Interval{N}, Interval{N}}}, ::AbstractHyperrectangle{N}) where {N}
convert(::Type{CartesianProductArray{N, Interval{N}}}, ::AbstractHyperrectangle{N}) where {N}
convert(::Type{MinkowskiSumArray}, ::MinkowskiSum{N, ST, MinkowskiSumArray{N, ST}}) where {N, ST}
convert(::Type{STAR}, ::AbstractPolyhedron{N}) where {N}
convert(::Type{STAR}, ::Star)
convert(::Type{SimpleSparsePolynomialZonotope}, ::AbstractZonotope)
convert(::Type{SimpleSparsePolynomialZonotope}, ::AbstractSparsePolynomialZonotope)
convert(::Type{SparsePolynomialZonotope}, ::AbstractZonotope{N}) where {N}
convert(::Type{SparsePolynomialZonotope}, ::SimpleSparsePolynomialZonotope{N}) where {N}
```
