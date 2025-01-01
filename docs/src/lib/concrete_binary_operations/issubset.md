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
w1 = ⊆(B1, B2, true)[2]
println(B1 ⊆ BI)
w2 = ⊆(B1, BI, true)[2]
println(B1 ⊆ H)
w3 = ⊆(B1, H, true)[2]
# 'B2 ⊆ B1' is not supported yet
# w11 = ⊆(B2, B1, true)[2]
println(B2 ⊆ BI)
w4 = ⊆(B2, BI, true)[2]
println(B2 ⊆ H)
println(BI ⊆ B1)
w5 = ⊆(BI, B1, true)[2]
println(BI ⊆ B2)
w6 = ⊆(BI, B2, true)[2]
println(BI ⊆ H)
w7 = ⊆(BI, H, true)[2]
println(H ⊆ B1)
w8 = ⊆(H, B1, true)[2]
println(H ⊆ B2)
w9 = ⊆(H, B2, true)[2]
println(H ⊆ BI)
w10 = ⊆(H, BI, true)[2];
```

```@example binary_set_operations
witnesses = [w1, w2, w3, w4, w5, w6, w7, w8, w9, w10]

plot1 = plot()
plot_sets(sets)
plot_points(witnesses, "w")
plot1
```

## Methods

```@docs
⊆(::LazySet, ::LazySet, ::Bool=false)
⊆(::LazySet, ::AbstractHyperrectangle, ::Bool=false)
⊆(::AbstractPolytope, ::AbstractHyperrectangle, ::Bool=false)
⊆(::AbstractZonotope, ::AbstractHyperrectangle)
⊆(::AbstractHyperrectangle, ::AbstractHyperrectangle, ::Bool=false)
⊆(::LazySet, ::AbstractPolyhedron, ::Bool=false)
⊆(::AbstractSingleton, ::LazySet, ::Bool=false)
⊆(::AbstractSingleton, ::AbstractHyperrectangle, ::Bool=false)
⊆(::AbstractSingleton, ::AbstractSingleton, ::Bool=false)
⊆(::Union{Ball2, Ballp}, ::AbstractSingleton, ::Bool=false)
⊆(::LineSegment, ::LazySet, ::Bool=false)
⊆(::LineSegment, ::AbstractHyperrectangle, ::Bool=false)
⊆(::Interval, ::UnionSet, ::Bool=false)
⊆(::LazySet, ::EmptySet, ::Bool=false)
⊆(::UnionSet, ::LazySet, ::Bool=false)
⊆(::UnionSetArray, ::LazySet, ::Bool=false)
⊆(::LazySet, ::Universe, ::Bool=false)
⊆(::Universe, ::LazySet, ::Bool=false)
⊆(::LazySet, ::Complement, ::Bool=false)
⊆(::CartesianProduct, ::CartesianProduct, ::Bool=false)
⊆(::CartesianProductArray, ::CartesianProductArray, ::Bool=false)
⊆(::AbstractZonotope, ::AbstractHyperrectangle, ::Bool=false)
⊆(::LazySet, ::UnionSetArray, ::Bool=false; ::Bool=true)
```

## Strict subset check

```@docs
⊂
```
