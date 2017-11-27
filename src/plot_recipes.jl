# ====================================
# Plot recipes for an abstract LazySet
# ====================================

"""
    plot_lazyset(X::T; ...) where {T<:LazySet}

Plot a lazy set in two-dimensions using an axis-aligned approximation.

### Input

- `X` -- a convex set

### Examples

```julia
julia> using LazySets, Plots
julia> B = BallInf(ones(2), 0.1)
julia> plot(2.0 * B)
```

### Notes

This recipe detects if the axis aligned approximation is such that the first two
vertices returned by `vertices_list` are the same. In that case, a scatter plot
is made (instead of the shape). This use case arises, for example, if you want
to plot operations on a singleton.
"""
@recipe function plot_lazyset(X::T; color="blue", label="",
                              grid=true, alpha=0.5) where {T<:LazySet}

    P = Approximations.overapproximate(X)
    vlist = hcat(vertices_list(P)...).'
    (x, y) = vlist[:, 1], vlist[:, 2]

    seriestype := norm(vlist[1, :] - vlist[2, :]) ≈ 0 ? :scatter : :shape

    x, y
end

"""
    plot_lazyset(X::Vector{T}) where {T<:LazySet}

Plot an array of lazy sets in two-dimensions using an axis-aligned approximation.

### Input

- `X` -- an array of convex sets

### Examples

```julia
julia> using LazySets, Plots
julia> B1 = BallInf(zeros(2), 0.4)
julia> B2 = BallInf(ones(2), 0.4)
julia> plot([B1, B2])
```
"""
@recipe function plot_lazyset(X::Vector{T}; seriescolor="blue", label="",
                              grid=true, alpha=0.5) where {T<:LazySet}

    seriestype := :shape

    for Xi in X
        Pi = Approximations.overapproximate(Xi)
        vlist = hcat(vertices_list(Pi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end

"""
    plot_lazyset(X::T, ε::Float64; ...) where {T<:LazySet}

Plot a lazy set in two-dimensions using iterative refinement.

### Input

- `X` -- a convex set
- `ε` -- approximation error bound

### Examples

```julia
julia> using LazySets, Plots
julia> B = BallInf(ones(2), 0.1)
julia> plot(randn(2, 2) * B, 1e-3)
```
"""
@recipe function plot_lazyset(X::T, ε::Float64; color="blue", label="",
                              grid=true, alpha=0.5) where {T<:LazySet}

    seriestype := :shape

    P = Approximations.overapproximate(X, ε)
    vlist = hcat(vertices_list(P)...).'
    (x, y) = vlist[:, 1], vlist[:, 2]

     x, y
end

"""
    plot_lazyset(X::Vector{T}, ε::Float64; ...) where {T<:LazySet}

Plot an array of lazy sets in two-dimensions using iterative refinement.

### Input

- `X` -- an array of convex sets
- `ε` -- approximation error bound

### Examples

```julia
julia> using LazySets, Plots
julia> B1 = BallInf(zeros(2), 0.4)
julia> B2 = Ball2(ones(2), 0.4)
julia> plot([B1, B2], 1e-4)
```
"""
@recipe function plot_lazyset(X::Vector{T}, ε::Float64; seriescolor="blue", label="",
                              grid=true, alpha=0.5) where {T<:LazySet}

    seriestype := :shape

    for Xi in X
        Pi = Approximations.overapproximate(Xi, ε)
        vlist = hcat(vertices_list(Pi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end

# =========================
# Plot recipes for polygons
# =========================

"""
    plot_polygon(P::Union{HPolygon, HPolygonOpt}; ...)

Plot a polygon given in constraint form.

### Input

- `P` -- a polygon in constraint representation

### Examples

```julia
julia> using LazySets, Plots
julia> P = HPolygon([LinearConstraint([1.0, 0.0], 0.6), LinearConstraint([0.0, 1.0], 0.6),
                     LinearConstraint([-1.0, 0.0], -0.4), LinearConstraint([0.0, -1.0], -0.4)])
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

Plot an array of polygons given in constraint form.

### Input

- `P` -- an array of polygons in constraint representation

### Examples

```julia
julia> using LazySets, Plots
julia> P1 = HPolygon([LinearConstraint([1.0, 0.0], 0.6), LinearConstraint([0.0, 1.0], 0.6),
                      LinearConstraint([-1.0, 0.0], -0.4), LinearConstraint([0.0, -1.0], -0.4)])
julia> P2 = HPolygon([LinearConstraint([2.0, 0.0], 0.6), LinearConstraint([0.0, 2.0], 0.6),
                      LinearConstraint([-2.0, 0.0], -0.4), LinearConstraint([0.0, -2.0], -0.4)])
julia> plot([P1, P2])
```
"""
@recipe function plot_polygons(P::Union{Vector{HPolygon}, Vector{HPolygonOpt}};
                              seriescolor="blue", label="", grid=true, alpha=0.5)

    seriestype := :shape

    for Pi in P
        vlist = hcat(vertices_list(Pi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end

"""
    plot_polygon(P::VPolygon; ...)

Plot a polygon given in vertex representation.

### Input

- `P` -- a polygon in vertex representation

### Examples

```julia
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

Plot an array of polygons given in vertex representation.

### Input

- `P` -- an array of polygons in vertex representation

### Examples

```julia
julia> using LazySets, Plots
julia> P1 = VPolygon([[0.6, 0.6], [0.4, 0.6], [0.4, 0.4], [0.6, 0.4]])
julia> P2 = VPolygon([[0.3, 0.3], [0.2, 0.3], [0.2, 0.2], [0.3, 0.2]])
julia> plot([P1, P2])
```
"""
@recipe function plot_polygons(P::Vector{VPolygon};
                              seriescolor="blue", label="", grid=true, alpha=0.5)

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

- `X` -- singleton, i.e. a one-element set

### Examples

```julia
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
    plot_singleton(X::Vector{Singleton}; ...)

Plot a list of singletons.

### Input

- `X` -- a list of singletons, i.e. a vector of one-element sets

### Examples

```julia
julia> using LazySets, Plots
julia> plot([Singleton([0.0, 0.0]), Singleton([1., 0]), Singleton([0.5, .5])])
```

Three-dimensional singletons can be plotted as well:

```julia
julia> using LazySets, Plots
julia> a, b, c = zeros(3), [1.0, 0, 0], [0.0, 1., 0];
julia> plot([Singleton(a), Singleton(b), Singleton(c)])
```
"""
@recipe function plot_singleton(X::Vector{Singleton};
                               color="blue", label="", grid=true,
                               legend=false)

    seriestype := :scatter

    [Tuple(Xi.element) for Xi in X]
end

# ============================
# Plot recipes for zonotopes
# ============================

"""
    plot_polygon(Z::Zonotope; ...)

Plot a zonotope, by enumerating its vertices.

### Input

- `Z` -- a zonotope represented by its center and generators

### Examples

```julia
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

- `Z` -- a linear array of zonotopes

### Examples

```julia
julia> using LazySets, Plots
julia> Z1 = Zonotope(zeros(2), [[0.6, 0.6], [0.4, 0.6], [0.4, 0.4], [0.6, 0.4]])
julia> Z2 = Zonotope(zeros(2), [[0.3, 0.3], [0.2, 0.3], [0.2, 0.2], [0.3, 0.2]])
julia> plot([Z1, Z2])
```
"""
@recipe function plot_zonotopes(V::Vector{Zonotope};
                              seriescolor="blue", label="", grid=true, alpha=0.5)

    seriestype := :shape

    for Zi in Z
        # we have to take the convex hull for the shape
        vlist = hcat(vertices_list(Zi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end
