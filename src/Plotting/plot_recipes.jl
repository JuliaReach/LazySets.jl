using RecipesBase
import RecipesBase.apply_recipe

using LazySets.Approximations: overapproximate, PolarDirections

# global values
DEFAULT_COLOR = :auto
DEFAULT_ALPHA = 0.5
DEFAULT_LABEL = ""
DEFAULT_GRID = true
DEFAULT_ASPECT_RATIO = :none
PLOT_PRECISION = 1e-3
PLOT_POLAR_DIRECTIONS = 40
DEFAULT_PLOT_LIMIT = 1000

function _extract_limits(p::RecipesBase.AbstractPlot)
    lims = Dict()
    if length(p) > 0
        subplot = p[1]
        for symbol in [:x, :y]
            lims[symbol] = subplot[Symbol(symbol,:axis)][:lims]
        end
    else
        lims[:x] = :auto
        lims[:y] = :auto
    end
    return lims
end

function _extract_extrema(p::RecipesBase.AbstractPlot)
    extrema = Dict()
    if length(p) > 0
        subplot = p[1]
        for symbol in [:x, :y]
            emin = subplot[Symbol(symbol,:axis)][:extrema].emin
            emax = subplot[Symbol(symbol,:axis)][:extrema].emax
            extrema[symbol] = (emin, emax)
        end
    else
        extrema[:x] = (-Inf, Inf)
        extrema[:y] = (-Inf, Inf)
    end
    return extrema
end

function _update_plot_limits!(lims, X::LazySet)
    box = box_approximation(X)
    box_min = low(box)
    box_max = high(box)
    for (idx,symbol) in enumerate([:x, :y])
        if lims[symbol] != :auto
            # scaling factor for obtaining offset from window width
            ϵ = 0.05
            # width of the plotting window in the direction `symbol`
            width = max(lims[symbol][2], box_max[idx]) -
                    min(lims[symbol][1], box_min[idx])

            # extend the current plot limits if the new set (plus a small offset) falls outside
            lims[symbol] = (min(lims[symbol][1], box_min[idx] - width*ϵ),
                            max(lims[symbol][2], box_max[idx] + width*ϵ))
        end
    end
    nothing
end

function _set_auto_limits_to_extrema!(lims, extr)
    # if the limit is :auto, set the limits to the current extrema
    for symbol in [:x, :y]
        if lims[symbol] == :auto
            emin, emax = extr[symbol]
            # if the extrema are (-Inf,Inf), i.e. an empty plot,
            # set it to (-1.,1.)
            lmin = isinf(emin) ? -1.0 : emin
            lmax = isinf(emax) ? 1.0 : emax
            lims[symbol] = (lmin, lmax)
        end
    end
    # otherwise keep the old limits
    nothing
end

function _bounding_hyperrectangle(lims, N)
    low_lim = [lims[:x][1] - DEFAULT_PLOT_LIMIT, lims[:y][1] - DEFAULT_PLOT_LIMIT]
    high_lim = [lims[:x][2] + DEFAULT_PLOT_LIMIT, lims[:y][2] + DEFAULT_PLOT_LIMIT]
    return Hyperrectangle(low=convert.(N,low_lim), high=convert.(N,high_lim))
end

"""
    plot_list(list::AbstractVector{VN}, [ε]::N=N(PLOT_PRECISION),
              [Nφ]::Int=PLOT_POLAR_DIRECTIONS, [fast]::Bool=false; ...)
        where {N<:Real, VN<:LazySet{N}}

Plot a list of convex sets.

### Input

- `list` -- list of convex sets (1D or 2D)
- `ε`    -- (optional, default: `PLOT_PRECISION`) approximation error bound
- `Nφ`   -- (optional, default: `PLOT_POLAR_DIRECTIONS`) number of polar
            directions (used to plot lazy intersections)
- `fast` -- (optional, default: `false`) switch for faster plotting but without
            individual plot recipes (see notes below)

### Notes

For each set in the list we apply an individual plot recipe.

The option `fast` provides access to a faster plotting scheme where all sets in
the list are first converted to polytopes and then plotted in one single run.
This, however, is not suitable when plotting flat sets (line segments,
singletons) because then the polytope plot recipe does not deliver good results.
Hence by default we do not use this option.
For plotting a large number of (non-flat) polytopes, we highly advise activating
this option.

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
@recipe function plot_list(list::AbstractVector{VN}, ε::N=N(PLOT_PRECISION),
                           Nφ::Int=PLOT_POLAR_DIRECTIONS, fast::Bool=false
                          ) where {N<:Real, VN<:LazySet{N}}
    if fast
        label --> DEFAULT_LABEL
        grid --> DEFAULT_GRID
        if DEFAULT_ASPECT_RATIO != :none
            aspect_ratio --> DEFAULT_ASPECT_RATIO
        end
        seriesalpha --> DEFAULT_ALPHA
        seriescolor --> DEFAULT_COLOR
        seriestype --> :shape

        first = true
        x = Vector{N}()
        y = Vector{N}()
        for Xi in list
            if Xi isa Intersection
                res = plot_recipe(Xi, ε, Nφ)
            else
                # hard-code overapproximation here to avoid individual
                # compilations for mixed sets
                Pi = overapproximate(Xi, ε)
                vlist = transpose(hcat(convex_hull(vertices_list(Pi))...))
                if isempty(vlist)
                    @warn "overapproximation during plotting was empty"
                    continue
                end
                res = vlist[:, 1], vlist[:, 2]
                # add first vertex to "close" the polygon
                push!(res[1], vlist[1, 1])
                push!(res[2], vlist[1, 2])
            end
            if isempty(res)
                continue
            else
                x_new, y_new = res
            end
            if first
                first = false
            else
                push!(x, N(NaN))
                push!(y, N(NaN))
            end
            append!(x, x_new)
            append!(y, y_new)
        end
        x, y
    else
        for Xi in list
            if Xi isa Intersection
                @series Xi, ε, Nφ
            else
                @series Xi, ε
            end
        end
    end
end

"""
    plot_lazyset(X::LazySet{N}, [ε]::N=N(PLOT_PRECISION); ...) where {N<:Real}

Plot a convex set.

### Input

- `X` -- convex set
- `ε` -- (optional, default: `PLOT_PRECISION`) approximation error bound

### Notes

See [`plot_recipe(::LazySet{<:Real})`](@ref).

For polyhedral set types (subtypes of `AbstractPolyhedron`), the argument `ε` is
ignored.

### Examples

```julia
julia> B = Ball2(ones(2), 0.1);

julia> plot(B, 1e-3)  # default accuracy value (explicitly given for clarity)

julia> plot(B, 1e-2)  # faster but less accurate than the previous call
```
"""
@recipe function plot_lazyset(X::LazySet{N}, ε::N=N(PLOT_PRECISION)
                             ) where {N<:Real}
    if dim(X) == 1
        plot_recipe(X, ε)
    else
        label --> DEFAULT_LABEL
        grid --> DEFAULT_GRID
        if DEFAULT_ASPECT_RATIO != :none
            aspect_ratio --> DEFAULT_ASPECT_RATIO
        end
        seriesalpha --> DEFAULT_ALPHA
        seriescolor --> DEFAULT_COLOR

        # extract limits and extrema of already plotted sets
        p = plotattributes[:plot_object]
        lims = _extract_limits(p)
        extr = _extract_extrema(p)

        if !isbounded(X)
            _set_auto_limits_to_extrema!(lims, extr)
            X = intersection(X, _bounding_hyperrectangle(lims, eltype(X)))

        # if there is already a plotted set and the limits are fixed,
        # automatically adjust the axis limits (e.g. after plotting a unbounded set)
        elseif length(p) > 0
            _update_plot_limits!(lims, X)
        end

        xlims --> lims[:x]
        ylims --> lims[:y]

        res = plot_recipe(X, ε)
        if isempty(res)
            res
        else
            x, y = res
            if length(x) == 1 || norm([x[1], y[1]] - [x[2], y[2]]) ≈ 0
                seriestype := :scatter
            else
                seriestype := :shape
            end
            x, y
        end
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
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR
    seriestype := :scatter

    # update manually set plot limits if necessary
    p = plotattributes[:plot_object]
    if length(p) > 0
        lims = _extract_limits(p)
        _update_plot_limits!(lims, S)
        xlims --> lims[:x]
        ylims --> lims[:y]
    end

    plot_recipe(S, ε)
end

"""
    plot_linesegment(X::Union{Interval{N}, LineSegment{N}}, [ε]::N=zero(N); ...)
        where {N<:Real}

Plot a line segment or an interval.

### Input

- `X` -- line segment or interval
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
@recipe function plot_linesegment(X::Union{Interval{N}, LineSegment{N}},
                                  ε::N=zero(N)) where {N<:Real}
    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    linecolor   --> DEFAULT_COLOR
    markercolor --> DEFAULT_COLOR
    markershape --> :circle
    seriestype := :path

    # update manually set plot limits if necessary
    p = plotattributes[:plot_object]
    if length(p) > 0
        lims = _extract_limits(p)
        _update_plot_limits!(lims, X)
        xlims --> lims[:x]
        ylims --> lims[:y]
    end

    plot_recipe(X, ε)
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
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end

    plot_recipe(∅)
end

"""
    plot_intersection(cap::Intersection{N}, [ε]::N=zero(N),
                      [Nφ]::Int=PLOT_POLAR_DIRECTIONS) where {N<:Real}

Plot a lazy intersection.

### Input

- `cap`  -- lazy intersection
- `ε`    -- (optional, default `0`) ignored, used for dispatch
- `Nφ`   -- (optional, default: `PLOT_POLAR_DIRECTIONS`) number of polar
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
                                   ε::N=zero(N),
                                   Nφ::Int=PLOT_POLAR_DIRECTIONS
                                  ) where {N<:Real}
    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR
    seriestype := :shape

    # extract limits and extrema of already plotted sets
    p = plotattributes[:plot_object]
    lims = _extract_limits(p)
    extr = _extract_extrema(p)

    if !isbounded(cap)
        _set_auto_limits_to_extrema!(lims, extr)
        bounding_box = _bounding_hyperrectangle(lims, eltype(cap))
        if !isbounded(cap.X)
            bounded_X = Intersection(bounding_box, cap.X)
            cap = Intersection(bounded_X, cap.Y)
        else
            cap = Intersection(bounding_box, cap.Y)
        end

    # if there is already a plotted set and the limits are fixed,
    # automatically adjust the axis limits (e.g. after plotting a unbounded set)
    elseif length(p) > 0
        _update_plot_limits!(lims, cap)
    end

    xlims --> lims[:x]
    ylims --> lims[:y]

    plot_recipe(cap, ε)
end
