using RecipesBase
import RecipesBase.apply_recipe

using LazySets.Approximations: overapproximate, PolarDirections

# global values
DEFAULT_COLOR = :auto
DEFAULT_ALPHA = 0.5
DEFAULT_LABEL = ""
DEFAULT_GRID = true
DEFAULT_ASPECT_RATIO = 1.0
DEFAULT_EPSILON = 1e-3
DEFAULT_POLAR_DIRECTIONS = 40

"""
    plot_list(list::AbstractVector{VN}, [ε]::N=N(DEFAULT_EPSILON),
              [Nφ]::Int=DEFAULT_POLAR_DIRECTIONS; ...) where {N<:Real,
                                                              VN<:LazySet{N}}

Plot a list of convex sets.

### Input

- `list` -- list of convex sets (1D or 2D)
- `ε`    -- (optional, default: `DEFAULT_EPSILON`) approximation error bound
- `Nφ`   -- (optional, default: `DEFAULT_POLAR_DIRECTIONS`) number of polar
            directions (used to plot lazy intersections)

### Notes

For each set in the list we apply an individual plot recipe.

### Examples

```julia
julia> B1 = BallInf(zeros(2), 0.4);

julia> B2 = BallInf(ones(2), 0.4);

julia> plot([B1, B2])
```

Some of the sets in the list may not be plotted precisely but rather
overapproximated first.
The second argument `ε` controls the accuracy of this overapproximation.

```julia
julia> Bs = [BallInf(zeros(2), 0.4), Ball2(ones(2), 0.4)];

julia> plot(Bs, 1e-3)  # default accuracy value (explicitly given for clarity)

julia> plot(Bs, 1e-2)  # faster but less accurate than the previous call
```
"""
@recipe function plot_list(list::AbstractVector{VN}, ε::N=N(DEFAULT_EPSILON),
                           Nφ::Int=DEFAULT_POLAR_DIRECTIONS
                          ) where {N<:Real, VN<:LazySet{N}}
    for Xi in list
        if Xi isa Intersection
            @series Xi, -one(N), Nφ
        else
            @series Xi, ε
        end
    end
end

"""
    plot_lazyset(X::LazySet{N}, [ε]::N=N(DEFAULT_EPSILON); ...) where {N<:Real}

Plot a convex set.

### Input

- `X` -- convex set
- `ε` -- (optional, default: `DEFAULT_EPSILON`) approximation error bound

### Notes

This recipe detects if the overapproximation contains duplicate vertices.
In that case, a scatter plot is used (instead of a shape plot).
This corner case arises, for example, from lazy linear maps of singletons.

Plotting of unbounded sets is not implemented yet (see
[#576](https://github.com/JuliaReach/LazySets.jl/issues/576)).

### Algorithm

We first assert that `X` is bounded.

One-dimensional sets are converted to an `Interval`.
Three-dimensional or higher-dimensional sets cannot be plotted.

For two-dimensional sets, we first compute a polygonal overapproximation.
The second argument, `ε`, corresponds to the error in Hausdorff distance between
the overapproximating set and `X`.
The default value `DEFAULT_EPSILON` is chosen such that the unit ball in the
2-norm is plotted with reasonable accuracy.
On the other hand, if you only want to produce a fast box-overapproximation of
`X`, pass `ε=Inf`.

In a second stage, we use the plot recipe for polygons.

### Examples

```julia
julia> B = Ball2(ones(2), 0.1);

julia> plot(B, 1e-3)  # default accuracy value (explicitly given for clarity)

julia> plot(B, 1e-2)  # faster but less accurate than the previous call
```
"""
@recipe function plot_lazyset(X::LazySet{N}, ε::N=N(DEFAULT_EPSILON)
                             ) where {N<:Real}
    @assert dim(X) <= 2 "cannot plot a $(dim(X))-dimensional set"
    @assert isbounded(X) "cannot plot an unbounded $(typeof(X))"

    if dim(X) == 1
        @series convert(Interval, X), ε
    else
        # construct epsilon-close polygon
        P = overapproximate(X, ε)
        # use polygon plot recipe
        @series P, ε
    end
end

"""
    plot_polyhedron(P::AbstractPolyhedron{N}, [ε]::N=zero(N); ...)
        where {N<:Real}

Plot a (bounded) polyhedron.

### Input

- `P` -- bounded polyhedron
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Algorithm

We first assert that `P` is bounded (i.e., that `P` is a polytope).

One-dimensional polytopes are converted to an `Interval`.
Three-dimensional or higher-dimensional polytopes are not supported.

For two-dimensional polytopes (i.e., polygons) we compute their set of vertices
using `vertices_list` and then plot the convex hull of these vertices.

### Examples

```julia
julia> P = HPolygon([LinearConstraint([1.0, 0.0], 0.6),
                     LinearConstraint([0.0, 1.0], 0.6),
                     LinearConstraint([-1.0, 0.0], -0.4),
                     LinearConstraint([0.0, -1.0], -0.4)]);

julia> plot(P)

julia> P = VPolygon([[0.6, 0.6], [0.4, 0.6], [0.4, 0.4], [0.6, 0.4]]);

julia> plot(P)
```
"""
@recipe function plot_polyhedron(P::AbstractPolyhedron{N}, ε::N=zero(N)
                                ) where {N<:Real}
    @assert dim(P) <= 2 "cannot plot a $(dim(P))-dimensional polyhedron"
    @assert isbounded(P) "cannot plot an unbounded $(typeof(P))"

    if dim(P) == 1
        @series convert(Interval, P), ε
    else
        label --> DEFAULT_LABEL
        grid --> DEFAULT_GRID
        aspect_ratio --> DEFAULT_ASPECT_RATIO
        seriesalpha --> DEFAULT_ALPHA
        seriescolor --> DEFAULT_COLOR

        vlist = transpose(hcat(convex_hull(vertices_list(P))...))
        if isempty(vlist)
            @warn "received a polyhedron with no vertices during plotting"
            return []
        end
        x, y = vlist[:, 1], vlist[:, 2]

        if length(x) == 1
            # a single point
            seriestype := :scatter
        else
            # add first vertex to "close" the polygon
            push!(x, vlist[1, 1])
            push!(y, vlist[1, 2])
            if norm(vlist[1, :] - vlist[2, :]) ≈ 0
                seriestype := :scatter
            else
                seriestype := :shape
            end
        end
        @series x, y
    end
end

"""
    plot_singleton(S::AbstractSingleton{N}, [ε]::N=zero(N); ...) where {N<:Real}

Plot a singleton.

### Input

- `S` -- singleton
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Examples

```julia
julia> plot(Singleton([0.5, 1.0]))
```
"""
@recipe function plot_singleton(S::AbstractSingleton{N}, ε::N=zero(N)
                               ) where {N<:Real}
    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    aspect_ratio --> DEFAULT_ASPECT_RATIO
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR
    seriestype := :scatter

    if dim(S) == 1
        @series [Tuple([element(S)[1], N(0)])]
    else
        @assert dim(S) ∈ [2, 3] "cannot plot a $(dim(S))-dimensional singleton"

        @series [Tuple(element(S))]
    end
end

"""
    plot_linesegment(L::LineSegment{N}, [ε]::N=zero(N); ...) where {N<:Real}

Plot a line segment.

### Input

- `L` -- line segment
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Examples

```julia
julia> L = LineSegment([0., 0.], [1., 1.]);

julia> plot(L)
```

To control the color of the line, use the `linecolor` keyword argument, and to
control the color of the end points, use the `markercolor` keyword argument.
To control the width, use `linewidth`.

```julia
julia> plot(L, markercolor="green", linecolor="red", linewidth=2.)
```

To omit the markers, use `markershape=:none`.
You also need to pass a value for `seriestype=:path` explicitly (this seems to
be an external bug).

```julia
julia> plot(L, seriestype=:path, markershape=:none)
```

A shorter alternative is to pass `marker=0`, but this may result in small dots
as markers based on the plotting backend.

```julia
julia> plot(L, marker=0)
```
"""
@recipe function plot_linesegment(L::LineSegment{N}, ε::N=zero(N)
                                 ) where {N<:Real}
    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    aspect_ratio --> DEFAULT_ASPECT_RATIO
    seriesalpha --> DEFAULT_ALPHA
    linecolor   --> DEFAULT_COLOR
    markercolor --> DEFAULT_COLOR
    markershape --> :circle
    seriestype := :path

    @series [Tuple(L.p); Tuple(L.q)]
end

"""
    plot_interval(I::Interval{N}, [ε]::N=zero(N); ...) where {N<:Real}

Plot an interval.

### Input

- `I` -- interval
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Notes

We convert the interval to a `LineSegment` with y coordinate equal to zero.
See the corresponding plot recipe for more discussion.

### Examples

```julia
julia> I = Interval(0.0, 1.0);

julia> plot(I)
```
"""
@recipe function plot_interval(I::Interval{N}, ε::N=zero(N)) where {N<:Real}
    @series LineSegment([min(I), N(0)], [max(I), N(0)]), ε
end

"""
    plot_emptyset(∅::EmptySet, [ε]::N=zero(N); ...)

Plot an empty set.

### Input

- `∅` -- empty set
- `ε` -- (optional, default: `0`) ignored, used for dispatch
"""
@recipe function plot_emptyset(∅::EmptySet{N}, ε::N=zero(N)) where {N<:Real}
    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    aspect_ratio --> DEFAULT_ASPECT_RATIO

    return []
end

"""
    plot_intersection(cap::Intersection{N}, [ε]::N=-one(N),
                      [Nφ]::Int=DEFAULT_POLAR_DIRECTIONS) where {N<:Real}

Plot a lazy intersection.

### Input

- `cap`  -- lazy intersection
- `ε`    -- (optional, default `-1`) ignored, used for dispatch
- `Nφ`   -- (optional, default: `DEFAULT_POLAR_DIRECTIONS`) number of polar
            directions used in the template overapproximation

### Notes

This function is separated from the main `LazySet` plot recipe because iterative
refinement is not available for lazy intersections (since it uses the support
vector (but see
[#1187](https://github.com/JuliaReach/LazySets.jl/issues/1187))).

Also note that if the set is a *nested* intersection, you may have to manually
overapproximate this set before plotting (see
`LazySets.Approximations.overapproximate` for details).

### Examples

```julia
julia> using LazySets.Approximations

julia> X = Ball2(zeros(2), 1.) ∩ Ball2(ones(2), 1.5);  # lazy intersection

julia> plot(X)
```

You can specify the accuracy of the overapproximation of the lazy intersection
by passing a higher value for `Nφ`, which stands for the number of polar
directions used in the overapproximation.
This number can also be passed to the `plot` function directly.

```julia
julia> plot(overapproximate(X, PolarDirections(100)))

julia> plot(X, -1., 100)  # equivalent to the above line
```
"""
@recipe function plot_intersection(cap::Intersection{N},
                                   ε::N=-one(N),
                                   Nφ::Int=DEFAULT_POLAR_DIRECTIONS
                                  ) where {N<:Real}
    if ε != -one(N)
        error("cannot plot a lazy intersection using iterative refinement " *
              "with accuracy threshold `ε = $ε`, because the exact support " *
              "vector of an intersection is currently not available; using " *
              "instead a set of `Nφ` template directions. To control the " *
              "number of directions, pass another integer argument as in " *
              "`plot(cap, -1., 40)`")
    end

    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    aspect_ratio --> DEFAULT_ASPECT_RATIO
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR

    if isempty(cap)
        return []
    elseif dim(cap) == 1
        @series convert(Interval, cap)
    else
        @assert dim(cap) == 2 "cannot plot a $(dim(S))-dimensional intersection"

        # construct polygon approximation using polar directions
        P = overapproximate(cap, PolarDirections{N}(Nφ))
        # use polygon plot recipe
        @series P, ε
    end
end
