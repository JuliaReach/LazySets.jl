```@contents
Pages = ["iterative_refinement.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets.Approximations
```

# Iterative refinement

```@docs
overapproximate_hausdorff
LocalApproximation
PolygonalOverapproximation
new_approx
addapproximation!
refine(::LocalApproximation, ::LazySet)
tohrep(::PolygonalOverapproximation)
convert(::Type{HalfSpace}, ::LocalApproximation)
```
