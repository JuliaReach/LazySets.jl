```@contents
Pages = ["issubset.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# Subset Check

The function `⊆` checks whether a set is a subset of another set.
It can optionally produce a witness if the subset relation does not hold.
The operator can be used in infix notation (`X ⊆ Y`).

!!! note
    `issubset` can be used as an alternative name to `⊆`.

## Examples

We use the following four sets for illustration.

```@example binary_set_operations
using LazySets, LazySets.Approximations, Plots
B1 = Ball1(-ones(2), 1.)
B2 = Ball2(ones(2), 1.)
BI = BallInf(zeros(2), 1.)
H = Hyperrectangle(ones(2), ones(2))
sets = [B1, B2, BI, H]

function plot_sets(sets)
    for S in sets
        println(S)
        plot!(S, 1e-2, fillalpha=0.1)
    end
end

function plot_points(points, prefix)
    for i in eachindex(points)
        p = points[i]
        num_occur = length(findfirst(x -> x == p, points[1:i]))
        x = p[1]
        y = p[2]
        if num_occur == 1
            x += 0.15
        elseif num_occur == 2
            y += 0.15
        elseif num_occur == 3
            x -= 0.15
        else
            y -= 0.15
        end
        plot!(Singleton(p))
        plot!(annotations=(x, y, text("$(prefix)$(i)")))
    end
end

plot1 = plot()
plot_sets(sets)
plot1
```

```@example binary_set_operations
println(B1 ⊆ B2)
w1 = issubset(B1, B2, true)[2]
println(B1 ⊆ BI)
w2 = issubset(B1, BI, true)[2]
println(B2 ⊆ B1)
w3 = issubset(B2, B1, true)[2]
println(B2 ⊆ BI)
w4 = issubset(B2, BI, true)[2]
println(B2 ⊆ H)
println(BI ⊆ B1)
w5 = issubset(BI, B1, true)[2]
println(BI ⊆ B2)
w6 = issubset(BI, B2, true)[2]
println(BI ⊆ H)
w7 = issubset(BI, H, true)[2]
println(H ⊆ B2)
w8 = issubset(H, B2, true)[2];
```

```@example binary_set_operations
witnesses = [w1, w2, w3, w4, w5, w6, w7, w8]

plot1 = plot(xlims=(-2, 2.3))
plot_sets(sets)
plot_points(witnesses, "w")
plot1
```

## Methods

```@meta
CurrentModule = LazySets.API
```

```@docs; canonical=false
issubset(::LazySet, ::LazySet)
```

```@meta
CurrentModule = LazySets
```

```@docs
issubset(::LazySet, ::AbstractHyperrectangle, ::Bool=false)
issubset(::AbstractPolytope, ::AbstractHyperrectangle, ::Bool=false)
issubset(::AbstractZonotope, ::AbstractHyperrectangle)
issubset(::LazySet, ::AbstractPolyhedron, ::Bool=false)
issubset(::AbstractSingleton, ::AbstractHyperrectangle, ::Bool=false)
issubset(::LineSegment, ::LazySet, ::Bool=false)
issubset(::LineSegment, ::AbstractHyperrectangle, ::Bool=false)
issubset(::Interval, ::UnionSet, ::Bool=false)
issubset(::LazySet, ::EmptySet, ::Bool=false)
issubset(::UnionSet, ::LazySet, ::Bool=false)
issubset(::UnionSetArray, ::LazySet, ::Bool=false)
issubset(::Universe, ::LazySet, ::Bool=false)
issubset(::LazySet, ::Complement, ::Bool=false)
issubset(::CartesianProduct, ::CartesianProduct, ::Bool=false)
issubset(::CartesianProductArray, ::CartesianProductArray, ::Bool=false)
issubset(::AbstractZonotope, ::AbstractHyperrectangle, ::Bool=false)
issubset(::LazySet, ::UnionSetArray, ::Bool=false; ::Bool=true)
```

## Strict subset check

```@meta
CurrentModule = LazySets.API
```

```@docs; canonical=false
⊂(::LazySet, ::LazySet)
```

```@meta
CurrentModule = LazySets
```

```@docs
⊂(::LazySet, ::LazySet)
```
