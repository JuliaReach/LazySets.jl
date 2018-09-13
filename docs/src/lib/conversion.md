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
convert(::Type{HPolygonOpt}, ::HPolygon)
convert(::Type{HPolytope}, ::HPolygon)
convert(::Type{HPolygon}, ::HPolytope)
convert(::Type{Zonotope}, ::AbstractHyperrectangle)
convert(::Type{Hyperrectangle}, ::Interval)
convert(::Type{HPolytope}, ::AbstractHyperrectangle)
convert(::Type{HPolygon}, ::AbstractHyperrectangle)
```
