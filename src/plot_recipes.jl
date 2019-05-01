using RecipesBase
import RecipesBase.apply_recipe

using LazySets.Approximations: overapproximate, PolarDirections

function warn_empty_polytope()
    @warn "received a polytope with no vertices during plotting"
end

# ====================================
# Plot recipes for an abstract LazySet
# ====================================

"""
    plot_lazyset(X::LazySet{N}, [ε]::N=N(1e-3); ...) where {N<:Real}

Plot a convex set in two dimensions.

### Input

- `X` -- convex set
- `ε` -- (optional, default: `1e-3`) approximation error bound

### Examples

```julia
julia> B = BallInf(ones(2), 0.1);

julia> plot(2.0 * B)
```

An iterative refinement method is applied to obtain an overapproximation of `X`
in constraint representation, which is then plotted. To improve the accuracy
of the iterative refinement, use the second argument using a small value:

```julia
julia> B = Ball2(ones(2), 0.1);

julia> plot(B, 1e-3);

julia> plot(B, 1e-2); # faster than the previous try, but less accurate
```

### Algorithm

In a first stage, an overapproximation of the given set to a polygon in constraint
representation is computed. The second argument, `ε`, corresponds to the error
in Hausdorff distance between the overapproximating set and `X`. The default
value `1e-3` is chosen such that the unit ball in the 2-norm is plotted with
reasonable accuracy. On the other hand, if you only want to produce a fast
box-overapproximation of `X`, pass `ε=Inf`.

In a second stage, the list of vertices of the overapproximation is computed
with the `vertices_list` function of the polygon.

A post-processing `convex_hull` is applied to the vertices list to ensure
that the shaded area inside the convex hull of the vertices is covered
correctly.

### Notes

This recipe detects if overapproximation is such that the first two vertices
returned by `vertices_list` are the same. In that case, a scatter plot is used
(instead of a shape plot). This corner case arises, for example, when lazy linear
maps of singletons.
"""
@recipe function plot_lazyset(X::LazySet{N}, ε::N=N(1e-3);
                              color="blue", label="", grid=true, alpha=0.5) where {N<:Real}

    @assert dim(X) == 2 "cannot plot a $(dim(X))-dimensional set using iterative refinement"

    P = overapproximate(X, ε)
    vlist = transpose(hcat(convex_hull(vertices_list(P))...))

    if isempty(vlist)
        warn_empty_polytope()
        return []
    end

    (x, y) = vlist[:, 1], vlist[:, 2]

    # add first vertex to "close" the polygon
    push!(x, vlist[1, 1])
    push!(y, vlist[1, 2])

    seriestype := norm(vlist[1, :] - vlist[2, :]) ≈ 0 ? :scatter : :shape

    x, y
end

"""
    plot_lazyset(X::Vector{XN}, ε::N=N(1e-3); ...) where {N<:Real, XN<:LazySet{N}}

Plot an array of convex sets in two dimensions.

### Input

- `X` -- array of convex sets
- `ε` -- (optional, default: `1e-3`) approximation error bound

### Examples

```julia
julia> B1 = BallInf(zeros(2), 0.4);

julia> B2 = BallInf(ones(2), 0.4);

julia> plot([B1, B2])
```

An iterative refinement method is applied to obtain an overapproximation of each
set in `X` in constraint representation, which is then plotted. To change the
default tolerance for the iterative refinement, use the second argument; it applies
to all sets in the array.

```julia
julia> B1 = BallInf(zeros(2), 0.4);

julia> B2 = Ball2(ones(2), 0.4);

julia> plot([B1, B2], 1e-1) # faster but less accurate 

julia> plot([B1, B2], 1e-4) # slower but more accurate 
```

### Algorithm

Refer to the documentation of `plot_lazyset(S::LazySet, [ε]::Float64; ...)` for
further details.
"""
@recipe function plot_lazysets(X::Vector{XN}, ε::N=N(1e-3);
                               seriescolor="blue", label="", grid=true,
                               alpha=0.5) where {N<:Real, XN<:LazySet{N}}

    seriestype := :shape

    for Xi in X
        if Xi isa EmptySet
            continue
        end
        @assert dim(Xi) == 2 "cannot plot a $(dim(Xi))-dimensional set using " *
                             "iterative refinement"
        Pi = overapproximate(Xi, ε)
        vlist = transpose(hcat(convex_hull(vertices_list(Pi))...))

        if isempty(vlist)
            warn_empty_polytope()
            continue
        end

        x, y = vlist[:, 1], vlist[:, 2]

        # add first vertex to "close" the polygon
        push!(x, vlist[1, 1])
        push!(y, vlist[1, 2])

        @series (x, y)
    end
end

# ==============================
# Plot recipes for 2D polytopes
# ==============================

"""
    plot_polytope(P::AbstractPolytope, [ε]::N=zero(N); ...) where {N<:Real}

Plot a 2D polytope as the convex hull of its vertices.

### Input

- `P` -- polygon or polytope
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Examples

```julia
julia> P = HPolygon([LinearConstraint([1.0, 0.0], 0.6),
                     LinearConstraint([0.0, 1.0], 0.6),
                     LinearConstraint([-1.0, 0.0], -0.4),
                     LinearConstraint([0.0, -1.0], -0.4)]);

julia> plot(P)
```

This recipe also applies if the polygon is given in vertex representation:
    
```julia
julia> P = VPolygon([[0.6, 0.6], [0.4, 0.6], [0.4, 0.4], [0.6, 0.4]]);

julia> plot(P);
```

### Algorithm

This function checks that the polytope received is two-dimensional, then
computes its vertices accessing its `vertices_list` function and takes their
convex hull. The set is plotted and shaded using the `:shape` series type from
`Plots`. 
"""
@recipe function plot_polytope(P::AbstractPolytope{N}, ε::N=zero(N);
                               color="blue", label="", grid=true, alpha=0.5) where {N<:Real}

    # for polytopes
    @assert dim(P) == 2 "cannot plot a $(dim(P))-dimensional polytope"
    seriestype := :shape

    points = convex_hull(vertices_list(P))
    vlist = transpose(hcat(points...))

    if isempty(vlist)
        warn_empty_polytope()
        return []
    end

    (x, y) = vlist[:, 1], vlist[:, 2]

    # add first vertex to "close" the polygon
    push!(x, vlist[1, 1])
    push!(y, vlist[1, 2])

    x, y
end

"""
    plot_polytopes(P::Vector{PN}, [ε]::N=zero(N); ...) where {N<:Real, PN<:AbstractPolytope{N}}

Plot an array of 2D polytopes as the convex hull of their vertices.

### Input

- `P` -- vector of polygons or polytopes
- `ε` -- (optional, default: `0.0`) ignored, used for dispatch

### Examples

```julia
julia> P = HPolygon([LinearConstraint([1.0, 0.0], 0.6),
                     LinearConstraint([0.0, 1.0], 0.6),
                     LinearConstraint([-1.0, 0.0], -0.4),
                     LinearConstraint([0.0, -1.0], -0.4)]);

julia> plot(P)
```

This recipe also applies if the polygon is given in vertex representation:
    
```julia
julia> P = VPolygon([[0.6, 0.6], [0.4, 0.6], [0.4, 0.4], [0.6, 0.4]]);

julia> plot(P)
```

### Algorithm

This function checks that the polytope received is two-dimensional, then
computes its vertices accessing its `vertices_list` function and takes their
convex hull. The set is plotted and shaded using the `:shape` series type from
`Plots`. 
"""
@recipe function plot_polytopes(P::Vector{PN}, ε::N=zero(N);
                                seriescolor="blue", label="", grid=true,
                                alpha=0.5) where {N<:Real, PN<:AbstractPolytope{N}}

    # it is assumed that the polytopes are two-dimensional
    seriestype := :shape

    for Pi in P
        @assert dim(Pi) == 2 "cannot plot a $(dim(Pi))-dimensional polytope"
        points = convex_hull(vertices_list(Pi))
        vlist = transpose(hcat(points...))

        if isempty(vlist)
            warn_empty_polytope()
            continue
        end

        x, y = vlist[:, 1], vlist[:, 2]

        # add first vertex to "close" the polygon
        push!(x, vlist[1, 1])
        push!(y, vlist[1, 2])

        @series (x, y)
    end
end

# ============================
# Plot recipes for singletons
# ============================

"""
    plot_singleton(S::AbstractSingleton{N}, [ε]::N=zero(N);; ...) where {N<:Real}

Plot a singleton.

### Input

- `S` -- singleton wrapping a point, i.e., a one-element set
- `ε` -- (optional, default: `0.0`) ignored, used for dispatch

### Examples

```julia
julia> plot(Singleton([0.5, 1.0]))
```
"""
@recipe function plot_singleton(S::AbstractSingleton{N}, ε::N=zero(N);
                                color="blue", label="", grid=true,
                                legend=false) where {N<:Real}

    seriestype := :scatter
    @assert dim(S) == 2 ||
            dim(S) == 3 "cannot plot a $(dim(S))-dimensional singleton"
    [Tuple(element(S))]
end

"""
    plot_singletons(S::Vector{SN}, ε::N=zero(N); ...) where {N<:Real, SN<:AbstractSingleton{N}}

Plot a list of singletons.

### Input

- `S` -- list of singletons, i.e., a vector of one-element sets
- `ε` -- (optional, default: `0.0`) ignored, used for dispatch

### Examples

```julia
julia> plot([Singleton([0.0, 0.0]), Singleton([1., 0]), Singleton([0.5, .5])])
```

Three-dimensional singletons can be plotted as well, with this same recipe:

```julia
julia> a, b, c = [0.0, 0, 0], [1.0, 0, 0], [0.0, 1., 0];

julia> plot([Singleton(a), Singleton(b), Singleton(c)])
```
"""
@recipe function plot_singletons(S::Vector{SN}, ε::N=zero(N);
                                 color="blue", label="",
                                 grid=true, legend=false) where {N<:Real, SN<:AbstractSingleton{N}}

    seriestype := :scatter

    if dim(S[1]) == 2
        @assert all([dim(p) == 2 for p in S]) "all points in this vector " *
            "should have the same dimension"
    elseif dim(S[1]) == 3
        @assert all([dim(p) == 3 for p in S]) "all points in this vector " *
            "should have the same dimension"
    else
        error("can only plot 2D or 3D vectors of singletons")
    end

    [Tuple(element(p)) for p in S]
end

# ==============================================
# Plot recipes for line segments and intervals
# ==============================================

"""
    plot_linesegment(L::LineSegment{N}, [ε]::N=zero(N); ...) where {N<:Real}

Plot a line segment.

### Input

- `L` -- line segment
- `ε` -- (optional, default: `0.0`) ignored, used for dispatch

### Examples

```julia
julia> L = LineSegment([0., 0.], [1., 1.]);

julia> plot(L)
```
"""
@recipe function plot_linesegment(L::LineSegment{N}, ε::N=zero(N); color="blue", label="",
                                  grid=true, alpha=0.5, legend=false,
                                  add_marker=true) where {N<:Real}

    seriestype := :path
    linecolor   --> color
    markershape --> (add_marker ? :circle : :none)
    markercolor --> color

    [Tuple(L.p); Tuple(L.q)]
end

"""
    plot_linesegments(L::Vector{LineSegment{N}}, [ε]::N=zero(N); ...) where {N<:Real}

Plot an array of line segments.

### Input

- `L` -- vector of line segments
- `ε` -- (optional, default: `0.0`) ignored, used for dispatch

### Examples

```julia
julia> L1 = LineSegment([0., 0.], [1., 1.]);

julia> L2 = LineSegment([1., 0.], [0., 1.]);

julia> plot([L1, L2])
```
"""
@recipe function plot_linesegments(L::Vector{LineSegment{N}}, ε::N=zero(N); color="blue",
                                   label="", grid=true, alpha=0.5, legend=false,
                                   add_marker=true) where {N<:Real}

    seriestype := :path
    linecolor   --> color
    markershape --> (add_marker ? :circle : :none)
    markercolor --> color

    for Li in L
        @series [Tuple(Li.p); Tuple(Li.q)]
    end
end

"""
    plot_interval(I::Interval{N}, [ε]::N=zero(N); ...) where {N<:Real}

Plot an interval.

### Input

- `I` -- interval
- `ε` -- (optional, default: `0.0`) ignored, used for dispatch

### Examples

```julia
julia> I = Interval(0.0, 1.0);

julia> plot(I)
```
"""
@recipe function plot_interval(I::Interval{N}, ε::N=zero(N); color=:auto, label="",
                               grid=true, alpha=0.5, legend=false, add_marker=true,
                               linewidth=2.) where {N<:Real}

    seriestype := :path
    linecolor   --> color
    markershape --> (add_marker ? :circle : :none)
    markercolor --> color

    [Tuple([min(I), 0.0]); Tuple([max(I), 0.0])]
end

"""
    plot_intervals(I::Vector{Interval{N}}, [ε]::N=zero(N); ...); ...) where {N<:Real}

Plot an array of intervals.

### Input

- `I` -- vector of intervals
- `ε` -- (optional, default: `0.0`) ignored, used for dispatch

### Examples

```julia
julia> I1 = Interval([0., 1.]);

julia> I2 = Interval([0.5, 2.]);

julia> plot([I1, I2])
```
"""
@recipe function plot_intervals(I::Vector{Interval{N}}, ε::N=zero(N); color=:auto, label="", grid=true,
                                alpha=0.5, legend=false, add_marker=true,
                                linewidth=2.0) where {N<:Real}

    seriestype := :path
    linecolor   --> color
    markershape --> (add_marker ? :circle : :none)
    markercolor --> color

    for Ii in I
        @series [Tuple([min(Ii), 0.0]); Tuple([max(Ii), 0.0])]
    end
end

# ==============================
# Plot recipe for the empty set
# ==============================

"""
    plot_emptyset(∅::EmptySet, [ε]; ...)

Plot an empty set.

### Input

- `∅` -- empty set
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

An empty figure.
"""
@recipe function plot_emptyset(∅::EmptySet{N}, ε::N=zero(N); label="", grid=true,
                               legend=false) where {N<:Real}
    return []
end

# ===============================
# Plot recipe for unbounded sets
# ===============================

"""
    plot_universe(U::Universe, [ε]; ...)

Plot the universal set.

### Input

- `U` -- universal set
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

An error; plotting the universal set is not implemented.
"""
@recipe function plot_universe(::Universe{N}, ε::N=zero(N)) where {N<:Real}
    error("cannot plot the universal set")
end

"""
    plot_line(L::Line, [ε]; ...)

Plot a line in 2D.

### Input

- `L` -- line
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

An error, since plotting a line is not implemented yet (see #576).
"""
@recipe function plot_line(::Line{N}, ε::N=zero(N)) where {N<:Real}
    error("cannot plot an infinite line")
end

"""
    plot_halfspace(H::HalfSpace, [ε]; ...)

Plot a half-space in 2D.

### Input

- `H` -- half-space
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

An error, since plotting a half-space is not implemented yet (see #576).
"""
@recipe function plot_halfspace(::HalfSpace{N}, ε::N=zero(N)) where {N<:Real}
    error("cannot plot a half-space")
end

"""
    plot_hyperplane(H::Hyperplane, [ε]; ...)

Plot a hyperplane in 2D.

### Input

- `H` -- hyperplane
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

An error, since plotting a hyperplane is not implemented yet (see #576).
"""
@recipe function plot_hyperplane(::Hyperplane{N}, ε::N=zero(N)) where {N<:Real}
    error("cannot plot a hyperplane")
end

# ===================================
# Plot recipe for lazy intersections
# ===================================

"""
    plot_intersection(X::Intersection{N}, [ε]::N=-one(N); Nφ=40,
                      color="blue", label="", grid=true, alpha=0.5) where {N<:Real}

Plot a lazy intersection.

### Input

- `X`  -- lazy intersection
- `ε`  -- (optional, default -1) ignored, used only for dispatch
- `Nφ` -- (optional, default: `40`) number of template directions used in the
          template overapproximation

### Output

A plot with the overapproximation of the given lazy intersection.

### Notes

This function is separated from the main one `LazySet` plot recipe because
the iterative refinement is not available for lazy intersections, since it uses
the support vector (but see #1187).

Also note that if the set is a *nested* intersection eg. the lazy linear map
a nested intersection you may have to manually overapproximate this set before
plotting, see `LazySets.Approximations.overapproximate` for details.

```julia
julia> using Polyhedra, LazySets.Approximations

julia> X = rand(Ball2) ∩ rand(Ball2); # lazy intersection

julia> plot(X)
```

The vertices list of an `HPolygon` has known issues (until #1306 is fixed).
Consider using the `Polyhedra` backend to compute the dual representation.
One can specify the accuaracy of the overapproximation of the lazy intersection
passing a higher value in `Nφ`, which stands for amount of polar directions chosen.

```julia
julia> Nφ = 100; # or a bigger number

julia> Po = overapproximate(X, PolarDirections(Nφ));

julia> P = convert(HPolytope, Po) # see issue #1306

julia> plot(P)
```
"""
@recipe function plot_intersection(X::Intersection{N}, ε::N=-one(N); Nφ=40,
                                   color="blue", label="", grid=true, alpha=0.5) where {N<:Real}

    @assert dim(X) == 2 "this recipe only plots two-dimensional sets"

    if ε != -one(N)
        error("cannot plot a lazy intersection using iterative refinement with " *
              "error threshold `ε = $ε`, because the exact support vector of an " *
              "intersection is not available; using instead a set of `Nφ` " *
              "template directions. To control the number of directions, pass the " *
              "variable Nφ as in `plot(X, Nφ=...)`")
    end

    P = convert(HPolygon, overapproximate(X, PolarDirections{N}(Nφ)))
    vlist = transpose(hcat(convex_hull(vertices_list(P))...))

    if isempty(vlist)
        warn_empty_polytope()
        return []
    end

    (x, y) = vlist[:, 1], vlist[:, 2]

    # add first vertex to "close" the polygon
    push!(x, vlist[1, 1])
    push!(y, vlist[1, 2])

    seriestype := norm(vlist[1, :] - vlist[2, :]) ≈ 0 ? :scatter : :shape

    x, y
end
