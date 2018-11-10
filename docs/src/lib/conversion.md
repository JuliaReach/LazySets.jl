# Conversion between set representations

This section of the manual lists the conversion functions between set
representations.

```@contents
Pages = ["conversion.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
DocTestSetup = quote
    using LazySets
end
```

```@docs
convert(T::Type{<:AbstractHPolygon}, P::VPolygon)
convert(::Type{Hyperrectangle}, ::Interval)
convert(::Type{HPolygon}, ::AbstractHyperrectangle)
convert(::Type{HPolygon}, ::HPolytope)
convert(::Type{HPolygonOpt}, ::HPolygon)
convert(::Type{HPolytope}, P::AbstractHPolygon)
convert(::Type{HPolytope}, ::AbstractHyperrectangle)
convert(::Type{HPolytope}, P::AbstractPolytope)
convert(::Type{HPolytope}, P::VPolytope)
convert(::Type{VPolygon}, P::AbstractHPolygon)
convert(::Type{VPolygon}, P::AbstractPolytope)
convert(::Type{VPolytope}, P::AbstractPolytope)
convert(::Type{VPolytope}, P::HPolytope)
convert(::Type{Zonotope}, ::AbstractHyperrectangle)
```
