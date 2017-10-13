# Common Set Representations

This section of the manual describes the basic set representation types.

```@contents
Pages = ["representations.md"]
```

```@meta
CurrentModule = LazySets
```

## Balls

Unit balls are defined by int center (vector) and radius (scalar), such as
infinity-norm balls,

$B_\infty(c, r) = \{ x ∈ \mathbb{R}^n : \Vert x - c\Vert_\infty \leq r \}.$

and Euclidean (2-norm) balls,

$B_2(c, r) = \{ x ∈ \mathbb{R}^n : \Vert x - c\Vert_2 \leq r \}.$

```@docs
Hyperrectangle
BallInf
Ball2
```

## Polygons

```@docs
HPolygon
HPolygonOpt
VPolygon
plot_polygon
tovrep
```

## More Types

```@docs
intersection
LinearConstraint
Line
Singleton
VoidSet
```