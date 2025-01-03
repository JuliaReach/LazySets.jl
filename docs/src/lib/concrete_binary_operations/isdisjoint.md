```@contents
Pages = ["isdisjoint.md"]
Depth = 3
```

```@meta
CurrentModule = LazySets
```

# Check for Disjointness of Sets

The function `isdisjoint` checks whether the intersection of two sets is empty.
It can optionally produce a witness if the intersection is nonempty.

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
println(isdisjoint(BI, H))
w1 = isdisjoint(BI, H, true)[2]
println(isdisjoint(B1, BI))
w2 = isdisjoint(B1, BI, true)[2]
println(isdisjoint(B1, H))
```

```@example binary_set_operations
witnesses = [w1, w2]

plot1 = plot()
plot_sets(sets)
plot_points(witnesses, "w")
plot1
```

## Methods

```@meta
CurrentModule = LazySets.API
```

```@docs; canonical=false
isdisjoint(::LazySet, ::LazySet)
```

```@meta
CurrentModule = LazySets
```

```@docs
isdisjoint(::LazySet, ::LazySet, ::Bool=false)
isdisjoint(::AbstractHyperrectangle, ::AbstractHyperrectangle, ::Bool=false)
isdisjoint(::LazySet, ::AbstractSingleton, ::Bool=false)
isdisjoint(::AbstractZonotope, ::Hyperplane, ::Bool=false)
isdisjoint(::LazySet, ::Hyperplane, ::Bool=false)
isdisjoint(::LazySet, ::HalfSpace, ::Bool=false)
isdisjoint(::AbstractPolyhedron, ::LazySet, ::Bool=false)
isdisjoint(::Complement, ::LazySet, ::Bool=false)
isdisjoint(::AbstractZonotope, ::AbstractZonotope, ::Bool=false)
isdisjoint(::CartesianProductArray, ::AbstractPolyhedron, ::Bool=false)
isdisjoint(::CartesianProductArray, ::CartesianProductArray, ::Bool=false)
isdisjoint(::CartesianProductArray, ::AbstractHyperrectangle, ::Bool=false)
```
