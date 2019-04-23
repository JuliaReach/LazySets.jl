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

```jldoctest test_plot_lazyset
julia> using Plots, LazySets

julia> B = BallInf(ones(2), 0.1);

julia> plot(2.0 * B);
```

An iterative refinement method is applied to obtain an overapproximation of `X`
in constraint representation, which is then plotted. To change the default tolerance
for the iterative refinement, use the second argument:

```jldoctest test_plot_lazyset
julia> B = Ball2(ones(2), 0.1);

julia> plot(B, 1e-3);

julia> plot(B, 1e-2); # faster than the previous algorithm, but less accurate
```

### Algorithm

In first stage, an overapproximation of the given set to a polygon in constraint
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
                              color="blue", label="", grid=true, alpha=0.5) where {N}

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
    plot_lazyset(Xk::Vector{LazySet{N}}, ε::N=N(1e-3); ...) where {N<:Real}

Plot an array of convex sets in two dimensions.

### Input

- `Xk` -- array of convex sets
- `ε` -- (optional, default: `1e-3`) approximation error bound

### Examples

```jldoctest
julia> using Plots, LazySets;

julia> B1 = BallInf(zeros(2), 0.4);

julia> B2 = BallInf(ones(2), 0.4);

julia> plot([B1, B2]);
```

An iterative refinement method is applied to obtain an overapproximation of each
set in `Xk` in constraint representation, which is then plotted. To change the
default tolerance for the iterative refinement, use the second argument:

```julia
```jldoctest
julia> using Plots, LazySets;

julia> B1 = BallInf(zeros(2), 0.4);

julia> B2 = Ball2(ones(2), 0.4);

julia> plot([B1, B2], 1e-1); # faster but less accurate 

julia> plot([B1, B2], 1e-4); # slower but more accurate 
```

### Algorithm

See the documentation of `plot_lazyset(S::LazySet, [ε]::Float64; ...)` for details.
"""
@recipe function plot_lazyset(Xk::Vector{LazySet{N}}, ε::N=N(1e-3);
                              seriescolor="blue", label="", grid=true,
                              alpha=0.5) where {N<:Real}

    seriestype := :shape

    for X in Xk
        if X isa EmptySet
            continue
        end
        @assert dim(X) == 2 "cannot plot a $(dim(X))-dimensional set using iterative refinement"
        Pi = overapproximate(X, ε)
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
    plot_polygon(P::AbstractPolytope; ...)

Plot a 2D polytope as the convex hull of its vertices.

### Input

- `P` -- polygon or polytope

### Examples

```jldoctest plotting_polytope
julia> using Plots, LazySets;

julia> P = HPolygon([LinearConstraint([1.0, 0.0], 0.6),
                     LinearConstraint([0.0, 1.0], 0.6),
                     LinearConstraint([-1.0, 0.0], -0.4),
                     LinearConstraint([0.0, -1.0], -0.4)]);

julia> plot(P);

```

This recipe also applies if the polygon is given in vertex representation:
    
```jldoctest plotting_polytope
julia> P = VPolygon([[0.6, 0.6], [0.4, 0.6], [0.4, 0.4], [0.6, 0.4]]);

julia> plot(P);

```
"""
@recipe function plot_polytope(P::AbstractPolytope;
                               color="blue", label="", grid=true, alpha=0.5)

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
    plot_polytopes(Xk::Vector{S}; ...)

Plot an array of 2D polytopes.

### Input

- `Xk` -- array of polytopes

### Examples

```jldoctest plotting_polytopes
julia> using Plots, LazySets;

julia> P1 = HPolygon([LinearConstraint([1.0, 0.0], 0.6),
                      LinearConstraint([0.0, 1.0], 0.6),
                      LinearConstraint([-1.0, 0.0], -0.4),
                      LinearConstraint([0.0, -1.0], -0.4)]);

julia> P2 = HPolygon([LinearConstraint([2.0, 0.0], 0.6),
                      LinearConstraint([0.0, 2.0], 0.6),
                      LinearConstraint([-2.0, 0.0], -0.4),
                      LinearConstraint([0.0, -2.0], -0.4)]);

julia> plot([P1, P2]);

```

```jldoctest plotting_polytopes
julia> P1 = VPolygon([[0.6, 0.6], [0.4, 0.6], [0.4, 0.4], [0.6, 0.4]]);

julia> P2 = VPolygon([[0.3, 0.3], [0.2, 0.3], [0.2, 0.2], [0.3, 0.2]]);

julia> plot([P1, P2]);

```

### Notes

It is assumed that the given vector of polytopes is two-dimensional.
"""
@recipe function plot_polytopes(Xk::Vector{S};
                               seriescolor="blue", label="", grid=true,
                               alpha=0.5) where {S<:AbstractPolytope}

    # it is assumed that the polytopes are two-dimensional
    seriestype := :shape

    for Pi in Xk
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
    plot_singleton(X::AbstractSingleton; ...)

Plot a singleton.

### Input

- `X` -- singleton, i.e., a one-element set

### Examples

```jldoctest
julia> using Plots, LazySets;

julia> plot(Singleton([0.5, 1.0]));

```
"""
@recipe function plot_singleton(point::AbstractSingleton;
                                color="blue", label="", grid=true,
                                legend=false)

    seriestype := :scatter
    @assert dim(point) == 2 ||
            dim(point) == 3 "cannot plot a $(dim(point))-dimensional singleton"
    [Tuple(element(point))]
end

"""
    plot_singleton(Xk::Vector{S}; ...) where {S<:AbstractSingleton}

Plot a list of singletons.

### Input

- `Xk` -- list of singletons, i.e., a vector of one-element sets

### Examples

```jldoctest
julia> using Plots, LazySets;

julia> plot([Singleton([0.0, 0.0]), Singleton([1., 0]), Singleton([0.5, .5])]);

```

Three-dimensional singletons can be plotted as well:

```jldoctest
julia> using Plots, LazySets;

julia> a, b, c = zeros(3), [1.0, 0, 0], [0.0, 1., 0];

julia> plot([Singleton(a), Singleton(b), Singleton(c)]);

```
"""
@recipe function plot_singleton(Xk::Vector{S};
                                color="blue", label="", grid=true, legend=false
                               ) where {S<:AbstractSingleton}

    seriestype := :scatter

    if dim(Xk[1]) == 2
        @assert all([dim(pi) == 2 for pi in Xk]) "all points in this vector " *
            "should have the same dimension"
    elseif dim(Xk[1]) == 3
        @assert all([dim(pi) == 3 for pi in Xk]) "all points in this vector " *
            "should have the same dimension"
    else
        error("can only plot 2D or 3D vectors of singletons")
    end

    [Tuple(element(point)) for point in Xk]
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
julia> using Plots, LazySets;

julia> L = LineSegment([0., 0.], [1., 1.]);

julia> plot(L);

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
    plot_linesegments(Xk::Vector{S}; ...) where {S<:LineSegment}

Plot an array of line segments.

### Input

- `Xk` -- linear array of line segments

### Examples

```jldoctest
julia> using Plots, LazySets;

julia> L1 = LineSegment([0., 0.], [1., 1.]);

julia> L2 = LineSegment([1., 0.], [0., 1.]);

julia> plot([L1, L2]);

```
"""
@recipe function plot_linesegments(Xk::Vector{S}; color="blue",
                                   label="", grid=true, alpha=0.5, legend=false,
                                   add_marker=true) where {S<:LineSegment}

    seriestype := :path
    linecolor   --> color
    markershape --> (add_marker ? :circle : :none)
    markercolor --> color

    for Li in Xk
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
julia> using Plots, LazySets;

julia> I = Interval(0.0, 1.0);

julia> plot(I);

```
"""
@recipe function plot_interval(I::Interval; color=:auto, label="", grid=true,
                               alpha=0.5, legend=false, add_marker=true,
                               linewidth=2.)

    seriestype := :path
    linecolor   --> color
    markershape --> (add_marker ? :circle : :none)
    markercolor --> color

    [Tuple([min(I), 0.0]); Tuple([max(I), 0.0])]
end

"""
    plot_intervals(Xk::Vector{S}; ...) where {S<:Interval}

Plot an array of intervals.

### Input

- `Xk` -- linear array of intervals

### Examples

```jldoctest
julia> using Plots, LazySets;

julia> I1 = Interval([0., 1.]);

julia> I2 = Interval([0.5, 2.]);

julia> plot([I1, I2]);

```
"""
@recipe function plot_intervals(Xk::Vector{S}; color=:auto, label="", grid=true,
                                alpha=0.5, legend=false, add_marker=true,
                                linewidth=2.0) where {S<:Interval}

    seriestype := :path
    linecolor   --> color
    markershape --> (add_marker ? :circle : :none)
    markercolor --> color

    for Ii in Xk
        @series [Tuple([min(Ii), 0.0]); Tuple([max(Ii), 0.0])]
    end
end

# ==============================
# Plot recipe for the empty set
# ==============================

"""
    plot_emptyset(∅::EmptySet, [ε::Float64=0.0]; ...)

Plot an empty set.

### Input

- `∅` -- empty set
- `ε` -- (optional, default: `0.0`) approximation error bound

### Examples

```jldoctest
julia> using Plots, LazySets;

julia> plot(∅);

julia> plot(∅, 1e-2);

```
"""
@recipe function plot_emptyset(∅::EmptySet{N}, ε::N=N(0.0); label="", grid=true,
                               legend=false) where {N<:Real}
    return []
end

@recipe function plot_universe(::Universe{N}, ε::N=N(0.0)) where {N<:Real}
    error("cannot plot the universal set")
end

@recipe function plot_line(::Line{N}, ε::N=N(0.0)) where {N<:Real}
    error("cannot plot an infinite line")
end

@recipe function plot_halfspace(::HalfSpace{N}, ε::N=N(0.0)) where {N<:Real}
    error("cannot plot a half-space")
end

@recipe function plot_hyperplane(::Hyperplane{N}, ε::N=N(0.0)) where {N<:Real}
    error("cannot plot a hyperplane")
end

@recipe function plot_intersection(X::Intersection; Nφ=10,
                                   color="blue", label="", grid=true, alpha=0.5)

    @assert dim(X) == 2 "this recipe only plots two-dimensional sets"

    P = convert(HPolygon, overapproximate(X, PolarDirections(Nφ)))
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

#@recipe function plot_intersection(::Intersection{N}, ε::N=N(0.0)) where {N<:Real}
#    error("cannot plot a lazy intersection using iterative refinement " *
#         "(the exact support vector of an intersection is not implemented)")
#end
