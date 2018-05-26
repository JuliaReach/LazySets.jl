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

### Algorithm

For any 2D lazy set we compute its box overapproximation, followed by the list of
vertices. A post-processing `convex_hull` is applied to the vertices list;
this ensures that the shaded area inside the convex hull of the vertices is covered
correctly.

### Notes

This recipe detects if the axis-aligned approximation is such that the first two
vertices returned by `vertices_list` are the same. In that case, a scatter plot
is used (instead of a shape plot). This use case arises, for example, when
plotting singletons.
"""
@recipe function plot_lazyset(S::LazySet;
                              color="blue", label="", grid=true, alpha=0.5)

    P = Approximations.overapproximate(S)
    vlist = copy(transpose(hcat(convex_hull(vertices_list(P))...)))
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

### Algorithm

For each 2D lazy set in the array we compute its box overapproximation, followed
by the list of vertices. A post-processing `convex_hull` is applied to the vertices list;
this ensures that the shaded area inside the convex hull of the vertices is covered
correctly.
"""
@recipe function plot_lazyset(arr::Vector{<:LazySet};
                              seriescolor="blue", label="", grid=true,
                              alpha=0.5)

    seriestype := :shape

    for S in arr
        Pi = Approximations.overapproximate(S)
        vlist = copy(transpose(hcat(convex_hull(vertices_list(Pi))...)))
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
    vlist = copy(transpose(hcat(vertices_list(P)...)))
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
        vlist = copy(transpose(hcat(vertices_list(Pi)...)))
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

    vlist = copy(transpose(hcat(vertices_list(P)...)))
    (x, y) = vlist[:, 1], vlist[:, 2]

     x, y
end

"""
    plot_polygons(P::Vector{<:AbstractHPolygon}; ...)

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
@recipe function plot_polygons(P::Vector{<:AbstractHPolygon};
                               seriescolor="blue", label="", grid=true,
                               alpha=0.5)

    seriestype := :shape

    for Pi in P
        vlist = copy(transpose(hcat(vertices_list(Pi)...)))
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

    vlist = copy(transpose(hcat(vertices_list(P)...)))
    (x, y) = vlist[:, 1], vlist[:, 2]

     x, y
end

"""
    plot_polygons(P::Vector{<:VPolygon}; ...)

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
@recipe function plot_polygons(P::Vector{<:VPolygon};
                               seriescolor="blue", label="", grid=true,
                               alpha=0.5)

    seriestype := :shape

    for Pi in P
        vlist = copy(transpose(hcat(vertices_list(Pi)...)))
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end

# ============================
# Plot recipes for singletons
# ============================

"""
    plot_singleton(X::AbstractSingleton; ...)

Plot a singleton.

### Input

- `X` -- singleton, i.e., a one-element set

### Examples

```jldoctest
julia> using LazySets, Plots
julia> plot(Singleton([0.5, 1.0]))
```
"""
@recipe function plot_singleton(point::AbstractSingleton;
                                color="blue", label="", grid=true,
                                legend=false)

    seriestype := :scatter

    [Tuple(element(point))]
end

"""
    plot_singleton(arr::Vector{<:AbstractSingleton}; ...)

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
@recipe function plot_singleton(arr::Vector{<:AbstractSingleton};
                                color="blue", label="", grid=true, legend=false)

    seriestype := :scatter

    [Tuple(element(point)) for point in arr]
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
    vlist = copy(transpose(hcat(vertices_list(Z)...)))
    (x, y) = vlist[:, 1], vlist[:, 2]

     x, y
end

"""
    plot_zonotopes(Z::Vector{<:Zonotope}; ...)

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
@recipe function plot_zonotopes(Z::Vector{<:Zonotope};
                                seriescolor="blue", label="", grid=true,
                                alpha=0.5)

    seriestype := :shape

    for Zi in Z
        # we have to take the convex hull for the shape
        vlist = copy(transpose(hcat(vertices_list(Zi)...)))
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end

# =====================================
# Plot recipes for lines and intervals
# =====================================

"""
    plot_linesegment(L::LineSegment; ...)

Plot a line segment.

### Input

- `L` -- line segment

### Examples

```jldoctest
julia> using LazySets, Plots
julia> L = LineSegment([0., 0.], [1., 1.])
julia> plot(L)
```
"""
@recipe function plot_linesegment(L::LineSegment; color="blue", label="",
                                  grid=true, alpha=0.5, legend=false,
                                  add_marker=true)

    seriestype := :path
    linecolor   --> color
    markershape --> (add_marker ? :circle : :none)
    markercolor --> color

    [Tuple(L.p); Tuple(L.q)]
end

"""
    plot_linesegments(L::Vector{<:LineSegment}; ...)

Plot an array of line segments.

### Input

- `L` -- linear array of line segments

### Examples

```jldoctest
julia> using LazySets, Plots
julia> L1 = LineSegment([0., 0.], [1., 1.])
julia> L2 = LineSegment([1., 0.], [0., 1.])
julia> plot([L1, L2])
```
"""
@recipe function plot_linesegments(L::Vector{<:LineSegment}; color="blue",
                                   label="", grid=true, alpha=0.5, legend=false,
                                   add_marker=true)

    seriestype := :path
    linecolor   --> color
    markershape --> (add_marker ? :circle : :none)
    markercolor --> color

    for Li in L
        @series [Tuple(Li.p); Tuple(Li.q)]
    end
end

"""
    plot_interval(I::Interval; ...)

Plot an interval.

### Input

- `I` -- interval

### Examples

```jldoctest
julia> using LazySets, Plots
julia> I = Interval(0.0, 1.0)
julia> plot(I)
```
"""
@recipe function plot_linesegment(I::Interval; color=:auto, label="",
                                  grid=true, alpha=0.5, legend=false,
                                  add_marker=true, linewidth=2.)

    seriestype := :path
    linecolor   --> color
    markershape --> (add_marker ? :circle : :none)
    markercolor --> color

    [Tuple([low(I), 0.0]); Tuple([high(I), 0.0])]
end

"""
    plot_intervals(I::Vector{<:Interval}; ...)

Plot an array of intervals.

### Input

- `I` -- linear array of intervals

### Examples

```jldoctest
julia> using LazySets, Plots
julia> I1 = Interval([0., 1.])
julia> I2 = Interval([0.5, 2.])
julia> plot(I1, I2])
```
"""
@recipe function plot_intervals(I::Vector{<:Interval}; color=:auto,
                                   label="", grid=true, alpha=0.5, legend=false,
                                   add_marker=true, linewidth=2.)

    seriestype := :path
    linecolor   --> color
    markershape --> (add_marker ? :circle : :none)
    markercolor --> color

    for Ii in I
        @series [Tuple([low(Ii), 0.0]); Tuple([high(Ii), 0.0])]
    end
end
