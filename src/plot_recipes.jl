# ====================================
# Plot recipes for an abstract LazySet
# ====================================

"""
    plot_lazyset(S::LazySet; ...)

Plot a convex set in two dimensions using an axis-aligned approximation.

### Input

- `S` -- convex set

### Examples

```jldoctest
julia> using LazySets, Plots
julia> B = BallInf(ones(2), 0.1)
julia> plot(2.0 * B)
```

### Notes

This recipe detects if the axis-aligned approximation is such that the first two
vertices returned by `vertices_list` are the same. In that case, a scatter plot
is used (instead of a shape plot). This use case arises, for example, when
plotting singletons.
"""
@recipe function plot_lazyset(S::LazySet;
                              color="blue", label="", grid=true, alpha=0.5)

    P = Approximations.overapproximate(S)
    vlist = hcat(vertices_list(P)...).'
    (x, y) = vlist[:, 1], vlist[:, 2]

    seriestype := norm(vlist[1, :] - vlist[2, :]) ≈ 0 ? :scatter : :shape

    x, y
end

"""
    plot_lazyset(arr::Vector{<:LazySet})

Plot an array of convex sets in two dimensions using an axis-aligned
approximation.

### Input

- `arr` -- array of convex sets

### Examples

```jldoctest
julia> using LazySets, Plots
julia> B1 = BallInf(zeros(2), 0.4)
julia> B2 = BallInf(ones(2), 0.4)
julia> plot([B1, B2])
```
"""
@recipe function plot_lazyset(arr::Vector{<:LazySet};
                              seriescolor="blue", label="", grid=true,
                              alpha=0.5)

    seriestype := :shape

    for S in arr
        Pi = Approximations.overapproximate(S)
        vlist = hcat(vertices_list(Pi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end

"""
    plot_lazyset(S::LazySet, ε::Float64; ...)

Plot a lazy set in two dimensions using iterative refinement.

### Input

- `S` -- convex set
- `ε` -- approximation error bound

### Examples

```jldoctest
julia> using LazySets, Plots
julia> B = BallInf(ones(2), 0.1)
julia> plot(randn(2, 2) * B, 1e-3)
```
"""
@recipe function plot_lazyset(S::LazySet, ε::Float64;
                              color="blue", label="", grid=true, alpha=0.5)

    seriestype := :shape

    P = Approximations.overapproximate(S, ε)
    vlist = hcat(vertices_list(P)...).'
    (x, y) = vlist[:, 1], vlist[:, 2]

    x, y
end

"""
    plot_lazyset(arr::Vector{<:LazySet}, ε::Float64; ...)

Plot an array of lazy sets in two dimensions using iterative refinement.

### Input

- `arr` -- array of convex sets
- `ε` -- approximation error bound

### Examples

```jldoctest
julia> using LazySets, Plots
julia> B1 = BallInf(zeros(2), 0.4)
julia> B2 = Ball2(ones(2), 0.4)
julia> plot([B1, B2], 1e-4)
```
"""
@recipe function plot_lazyset(arr::Vector{<:LazySet}, ε::Float64;
                              seriescolor="blue", label="", grid=true,
                              alpha=0.5)

    seriestype := :shape

    for S in arr
        Pi = Approximations.overapproximate(S, ε)
        vlist = hcat(vertices_list(Pi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end

# =========================
# Plot recipes for polygons
# =========================

"""
    plot_polygon(P::Union{HPolygon, HPolygonOpt}; ...)

Plot a polygon in constraint representation.

### Input

- `P` -- polygon in constraint representation

### Examples

```jldoctest
julia> using LazySets, Plots
julia> P = HPolygon([LinearConstraint([1.0, 0.0], 0.6),
                     LinearConstraint([0.0, 1.0], 0.6),
                     LinearConstraint([-1.0, 0.0], -0.4),
                     LinearConstraint([0.0, -1.0], -0.4)])
julia> plot(P)
```
"""
@recipe function plot_polygon(P::Union{HPolygon, HPolygonOpt};
                              color="blue", label="", grid=true, alpha=0.5)

    seriestype := :shape

    vlist = hcat(vertices_list(P)...).'
    (x, y) = vlist[:, 1], vlist[:, 2]

     x, y
end

"""
    plot_polygons(P::Union{Vector{HPolygon}, Vector{HPolygonOpt}}; ...)

Plot an array of polygons in constraint representation.

### Input

- `P` -- array of polygons in constraint representation

### Examples

```jldoctest
julia> using LazySets, Plots
julia> P1 = HPolygon([LinearConstraint([1.0, 0.0], 0.6),
                      LinearConstraint([0.0, 1.0], 0.6),
                      LinearConstraint([-1.0, 0.0], -0.4),
                      LinearConstraint([0.0, -1.0], -0.4)])
julia> P2 = HPolygon([LinearConstraint([2.0, 0.0], 0.6),
                      LinearConstraint([0.0, 2.0], 0.6),
                      LinearConstraint([-2.0, 0.0], -0.4),
                      LinearConstraint([0.0, -2.0], -0.4)])
julia> plot([P1, P2])
```
"""
@recipe function plot_polygons(P::Union{Vector{HPolygon}, Vector{HPolygonOpt}};
                               seriescolor="blue", label="", grid=true,
                               alpha=0.5)

    seriestype := :shape

    for Pi in P
        vlist = hcat(vertices_list(Pi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end

"""
    plot_polygon(P::VPolygon; ...)

Plot a polygon in vertex representation.

### Input

- `P` -- polygon in vertex representation

### Examples

```jldoctest
julia> using LazySets, Plots
julia> P = VPolygon([[0.6, 0.6], [0.4, 0.6], [0.4, 0.4], [0.6, 0.4]])
julia> plot(P)
```
"""
@recipe function plot_polygon(P::VPolygon;
                              color="blue", label="", grid=true, alpha=0.5)

    seriestype := :shape

    vlist = hcat(vertices_list(P)...).'
    (x, y) = vlist[:, 1], vlist[:, 2]

     x, y
end

"""
    plot_polygons(P::Vector{VPolygon}; ...)

Plot an array of polygons in vertex representation.

### Input

- `P` -- array of polygons in vertex representation

### Examples

```jldoctest
julia> using LazySets, Plots
julia> P1 = VPolygon([[0.6, 0.6], [0.4, 0.6], [0.4, 0.4], [0.6, 0.4]])
julia> P2 = VPolygon([[0.3, 0.3], [0.2, 0.3], [0.2, 0.2], [0.3, 0.2]])
julia> plot([P1, P2])
```
"""
@recipe function plot_polygons(P::Vector{VPolygon};
                               seriescolor="blue", label="", grid=true,
                               alpha=0.5)

    seriestype := :shape

    for Pi in P
        vlist = hcat(vertices_list(Pi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end

# ============================
# Plot recipes for singletons
# ============================

"""
    plot_singleton(X::Singleton; ...)

Plot a singleton.

### Input

- `X` -- singleton, i.e., a one-element set

### Examples

```jldoctest
julia> using LazySets, Plots
julia> plot(Singleton([0.5, 1.0]))
```
"""
@recipe function plot_singleton(X::Singleton;
                                color="blue", label="", grid=true,
                                legend=false)

    seriestype := :scatter

    [Tuple(X.element)]
end

"""
    plot_singleton(arr::Vector{Singleton}; ...)

Plot a list of singletons.

### Input

- `arr` -- list of singletons, i.e., a vector of one-element sets

### Examples

```jldoctest
julia> using LazySets, Plots
julia> plot([Singleton([0.0, 0.0]), Singleton([1., 0]), Singleton([0.5, .5])])
```

Three-dimensional singletons can be plotted as well:

```jldoctest
julia> using LazySets, Plots
julia> a, b, c = zeros(3), [1.0, 0, 0], [0.0, 1., 0];
julia> plot([Singleton(a), Singleton(b), Singleton(c)])
```
"""
@recipe function plot_singleton(arr::Vector{Singleton};
                                color="blue", label="", grid=true, legend=false)

    seriestype := :scatter

    [Tuple(S.element) for S in arr]
end

# ============================
# Plot recipes for zonotopes
# ============================

"""
    plot_polygon(Z::Zonotope; ...)

Plot a zonotope by enumerating its vertices.

### Input

- `Z` -- zonotope

### Examples

```jldoctest
julia> using LazySets, Plots
julia> Z = Zonotope(ones(2), 0.2*[[1., 0], [0., 1], [1, 1]])
julia> plot(Z)
```
"""
@recipe function plot_zonotope(Z::Zonotope;
                               color="blue", label="", grid=true, alpha=0.5)

    seriestype := :shape

    # we have to take the convex hull for the shape
    vlist = hcat(vertices_list(Z)...).'
    (x, y) = vlist[:, 1], vlist[:, 2]

     x, y
end

"""
    plot_zonotopes(Z::Vector{Zonotope}; ...)

Plot an array of zonotopes.

### Input

- `Z` -- linear array of zonotopes

### Examples

```jldoctest
julia> using LazySets, Plots
julia> Z1 = Zonotope(zeros(2), [[0.6, 0.6], [0.4, 0.6], [0.4, 0.4], [0.6, 0.4]])
julia> Z2 = Zonotope(zeros(2), [[0.3, 0.3], [0.2, 0.3], [0.2, 0.2], [0.3, 0.2]])
julia> plot([Z1, Z2])
```
"""
@recipe function plot_zonotopes(V::Vector{Zonotope};
                                seriescolor="blue", label="", grid=true,
                                alpha=0.5)

    seriestype := :shape

    for Zi in Z
        # we have to take the convex hull for the shape
        vlist = hcat(vertices_list(Zi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end
