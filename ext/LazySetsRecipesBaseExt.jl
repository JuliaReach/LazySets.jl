module LazySetsRecipesBaseExt

using LazySets
using LazySets: plot_vlist, _plot_recipe_3d_polytope
using ReachabilityBase.Comparison: isapproxzero, _ztol
using RecipesBase: AbstractPlot, @recipe, @series
import LazySets: plot_recipe
import RecipesBase: apply_recipe  # required for Documenter to find docstrings

# global values
DEFAULT_COLOR = :auto
DEFAULT_ALPHA = 0.5
DEFAULT_LABEL = ""
DEFAULT_GRID = true
DEFAULT_ASPECT_RATIO = :none
PLOT_PRECISION = 1e-3
PLOT_POLAR_DIRECTIONS = 40
DEFAULT_PLOT_LIMIT = 1000

function _extract_limits(p::AbstractPlot, plotattributes::AbstractDict)
    lims = Dict()
    if length(p) > 0
        subplot = p[1]
        for symbol in [:x, :y]
            lims[symbol] = subplot[Symbol(symbol, :axis)][:lims]
        end
    else
        lims[:x] = :auto
        lims[:y] = :auto
    end

    # check whether the current call to `plot`/`plot!` passed new bounds
    if haskey(plotattributes, :xlims)
        lims[:x] = plotattributes[:xlims]
    end
    if haskey(plotattributes, :ylims)
        lims[:y] = plotattributes[:ylims]
    end

    return lims
end

function _extract_extrema(p::AbstractPlot)
    extrema = Dict()
    if length(p) > 0
        subplot = p[1]
        for symbol in [:x, :y]
            bounds = subplot[Symbol(symbol, :axis)][:extrema]
            emin = bounds.emin
            emax = bounds.emax
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
    isempty(box) && return nothing  # can happen if X is empty or flat

    box_min = low(box)
    box_max = high(box)
    n = dim(box)
    for (idx, symbol) in enumerate((:x, :y, :z))
        if idx > n
            break
        end
        if lims[symbol] != :auto
            # width of the plotting window including the new set in the `symbol` direction
            dmin = min(lims[symbol][1], box_min[idx])
            dmax = max(lims[symbol][2], box_max[idx])
            width = dmax - dmin

            # scaling factor for beautification
            ϵ = 0.03
            offset = width * ϵ

            # extend the current plot limits if the new set (plus a small offset) falls outside
            lims[symbol] = (min(lims[symbol][1], box_min[idx] - offset),
                            max(lims[symbol][2], box_max[idx] + offset))
        end
    end
    return nothing
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
    return nothing
end

function _bounding_hyperrectangle(lims, n, N)
    if n < 1 || n > 3
        throw(ArgumentError("cannot plot a $n-dimensional set"))
    end

    low_lim = Vector{N}(undef, n)
    high_lim = Vector{N}(undef, n)
    @inbounds for (i, s) in enumerate((:x, :y, :z))
        if i > n
            break
        end
        low_lim[i] = lims[s][1] - DEFAULT_PLOT_LIMIT
        high_lim[i] = lims[s][2] + DEFAULT_PLOT_LIMIT
    end
    return Hyperrectangle(; low=low_lim, high=high_lim)
end

"""
    plot_list(list::AbstractVector{VN}, [ε]::Real=N(PLOT_PRECISION),
              [Nφ]::Int=PLOT_POLAR_DIRECTIONS; [same_recipe]=false; ...)
        where {N, VN<:LazySet{N}}

Plot a list of sets.

### Input

- `list` -- list of sets (1D or 2D)
- `ε`    -- (optional, default: `PLOT_PRECISION`) approximation error bound
- `Nφ`   -- (optional, default: `PLOT_POLAR_DIRECTIONS`) number of polar
            directions (used to plot lazy intersections)
- `same_recipe` -- (optional, default: `false`) switch for faster plotting but
            without individual plot recipes (see notes below)

### Notes

For each set in the list we apply an individual plot recipe.

The option `same_recipe` provides access to a faster plotting scheme where all
sets in the list are first converted to polytopes and then plotted in one single
run. This, however, is not suitable when plotting flat sets (line segments,
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
@recipe function plot_list(list::AbstractVector{VN}, ε::Real=N(PLOT_PRECISION),  # COV_EXCL_LINE
                           Nφ::Int=PLOT_POLAR_DIRECTIONS;
                           same_recipe=false) where {N,VN<:LazySet{N}}
    if same_recipe
        label --> DEFAULT_LABEL
        grid --> DEFAULT_GRID
        if DEFAULT_ASPECT_RATIO != :none
            aspect_ratio --> DEFAULT_ASPECT_RATIO
        end
        seriesalpha --> DEFAULT_ALPHA
        seriescolor --> DEFAULT_COLOR
        seriestype --> :shape
        return _plot_list_same_recipe(list, ε, Nφ)
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

function _plot_list_same_recipe(list::AbstractVector{VN}, ε::Real=N(PLOT_PRECISION),
                                Nφ::Int=PLOT_POLAR_DIRECTIONS) where {N,VN<:LazySet{N}}
    first = true
    x = Vector{N}()
    y = Vector{N}()
    for Xi in list
        if Xi isa Intersection
            res = plot_recipe(Xi, ε, Nφ)
        elseif dim(Xi) == 1
            res = plot_recipe(Xi, ε)
        else
            # hard-code overapproximation here to avoid individual
            # compilations for mixed sets
            vlist = plot_vlist(Xi, ε)
            if isempty(vlist)
                @warn "overapproximation during plotting was empty"
                continue
            end
            vlist = transpose(hcat(vlist...))  # transpose vertices
            res = vlist[:, 1], vlist[:, 2]
            if length(res[1]) > 2
                # add first vertex to "close" the polygon
                push!(res[1], vlist[1, 1])
                push!(res[2], vlist[1, 2])
            end
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
    return x, y
end

# recipe for vector of singletons
@recipe function plot_list(list::AbstractVector{SN}) where {N,SN<:AbstractSingleton{N}}  # COV_EXCL_LINE
    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR
    seriestype --> :scatter

    return _plot_singleton_list(list)
end

# plot recipe for the union of singletons
@recipe function plot_list(X::UnionSetArray{N,SN}) where {N,SN<:AbstractSingleton{N}}  # COV_EXCL_LINE
    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR
    seriestype --> :scatter

    list = array(X)
    return _plot_singleton_list(list)
end

function _plot_singleton_list(list)
    n = dim(first(list))
    if n == 1
        _plot_singleton_list_1D(list)
    elseif n == 2
        _plot_singleton_list_2D(list)
    else
        throw(ArgumentError("plotting singletons is only available for " *
                            "dimensions one or two, but got dimension $n"))
    end
end

function _plot_singleton_list_1D(list::AbstractVector{SN}) where {N,SN<:AbstractSingleton{N}}
    m = length(list)

    x = Vector{N}(undef, m)
    y = zeros(N, m)

    @inbounds for (i, Xi) in enumerate(list)
        p = center(Xi)
        x[i] = p[1]
    end
    return x, y
end

function _plot_singleton_list_2D(list::AbstractVector{SN}) where {N,SN<:AbstractSingleton{N}}
    m = length(list)
    x = Vector{N}(undef, m)
    y = Vector{N}(undef, m)

    @inbounds for (i, Xi) in enumerate(list)
        p = center(Xi)
        x[i] = p[1]
        y[i] = p[2]
    end
    return x, y
end

"""
    plot_lazyset(X::LazySet{N}, [ε]::Real=N(PLOT_PRECISION); ...) where {N}

Plot a set.

### Input

- `X` -- set
- `ε` -- (optional, default: `PLOT_PRECISION`) approximation error bound

### Notes

This recipe just defines the default plotting options and then calls the
function `plot_recipe`, which then implements the set-specific plotting.

The argument `ε` is ignored by some set types, e.g., for polyhedra (subtypes of
`AbstractPolyhedron`).

### Examples

```julia
julia> B = Ball2(ones(2), 0.1);

julia> plot(B, 1e-3)  # default accuracy value (explicitly given for clarity here)

julia> plot(B, 1e-2)  # faster but less accurate than the previous call
```
"""
@recipe function plot_lazyset(X::LazySet{N}, ε::Real=N(PLOT_PRECISION)) where {N}  # COV_EXCL_LINE
    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end
    seriesalpha --> DEFAULT_ALPHA
    seriescolor --> DEFAULT_COLOR

    n = dim(X)
    if n == 1
        if !isbounded(X)
            # extract limits and extrema of already plotted sets
            p = plotattributes[:plot_object]
            lims = _extract_limits(p, plotattributes)
            extr = _extract_extrema(p)
            _set_auto_limits_to_extrema!(lims, extr)
            X = intersection(X, _bounding_hyperrectangle(lims, n, eltype(X)))
        end
        x, y = plot_recipe(X, ε)
        if length(x) == 1
            seriestype := :scatter
        end
        x, y
    elseif n == 2
        # extract limits and extrema of already plotted sets
        p = plotattributes[:plot_object]
        lims = _extract_limits(p, plotattributes)
        extr = _extract_extrema(p)

        if !isbounded(X)
            _set_auto_limits_to_extrema!(lims, extr)
            X = intersection(X, _bounding_hyperrectangle(lims, n, eltype(X)))

            # if there is already a plotted set and the limits are fixed,
            # automatically adjust the axis limits (e.g. after plotting an
            # unbounded set)
        elseif !isempty(X)
            _update_plot_limits!(lims, X)
        end

        xlims --> lims[:x]
        ylims --> lims[:y]

        res = plot_recipe(X, ε)
        if isempty(res)
            res
        else
            x, y = res
            if length(x) == 1 ||
               (length(x) == 2 && isapproxzero(norm([x[1], y[1]] - [x[2], y[2]])))
                # single point
                seriestype := :scatter
            elseif length(x) == 2
                # flat line segment
                linecolor --> DEFAULT_COLOR
                seriestype := :path
            else
                seriestype := :shape
            end
            x, y
        end
    elseif n == 3
        res = plot_recipe(X, ε)
        if isempty(res)
            res
        else
            x, y, z, i, j, k = res
            seriestype := :mesh3d
            connections := (i, j, k)
            x, y, z
        end
    else
        throw(ArgumentError("cannot plot a $n-dimensional set"))
    end
end

"""
    plot_singleton(S::AbstractSingleton{N}, [ε]::Real=zero(N); ...) where {N}

Plot a singleton.

### Input

- `S` -- singleton
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Examples

```julia
julia> plot(Singleton([0.5, 1.0]))
```
"""
@recipe function plot_singleton(S::AbstractSingleton{N}, ε::Real=zero(N)) where {N}  # COV_EXCL_LINE
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
        lims = _extract_limits(p, plotattributes)
        _update_plot_limits!(lims, S)
        xlims --> lims[:x]
        ylims --> lims[:y]
    end

    return plot_recipe(S, ε)
end

@recipe function plot_emptyset(∅::EmptySet{N}, ::Real=zero(N)) where {N}  # COV_EXCL_LINE
    label --> DEFAULT_LABEL
    grid --> DEFAULT_GRID
    if DEFAULT_ASPECT_RATIO != :none
        aspect_ratio --> DEFAULT_ASPECT_RATIO
    end

    return plot_recipe(∅)
end

"""
    plot_intersection(cap::Intersection{N}, [ε]::Real=zero(N),
                      [Nφ]::Int=PLOT_POLAR_DIRECTIONS) where {N}

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
overapproximate this set before plotting (see `overapproximate` for details).

### Examples

```julia
julia> X = Ball2(zeros(2), 1.) ∩ Ball2(ones(2), 1.5);  # lazy intersection

julia> plot(X)
```

You can specify the accuracy of the overapproximation of the lazy intersection
by passing an explicit value for `Nφ`, which stands for the number of polar
directions used in the overapproximation.
This number can also be passed to the `plot` function directly.

```julia
julia> plot(overapproximate(X, PolarDirections(100)))

julia> plot(X, 0.0, 100)  # equivalent to the above line
```
"""
@recipe function plot_intersection(cap::Intersection{N}, ε::Real=zero(N),  # COV_EXCL_LINE
                                   ::Int=PLOT_POLAR_DIRECTIONS) where {N}
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
    lims = _extract_limits(p, plotattributes)
    extr = _extract_extrema(p)

    n = dim(cap)
    if !isbounded(cap)
        _set_auto_limits_to_extrema!(lims, extr)
        bounding_box = _bounding_hyperrectangle(lims, n, eltype(cap))
        if !isbounded(first(cap))
            bounded_X = Intersection(bounding_box, first(cap))
            cap = Intersection(bounded_X, second(cap))
        else
            cap = Intersection(bounding_box, second(cap))
        end

        # if there is already a plotted set and the limits are fixed,
        # automatically adjust the axis limits (e.g. after plotting a unbounded set)
    elseif length(p) > 0
        _update_plot_limits!(lims, cap)
    end

    xlims --> lims[:x]
    if n > 1
        ylims --> lims[:y]
    end

    return plot_recipe(cap, ε)
end

# non-convex sets

@recipe function plot_union(cup::Union{UnionSet{N},UnionSetArray{N}}, ε::Real=N(PLOT_PRECISION);  # COV_EXCL_LINE
                            same_recipe=false, Nφ=PLOT_POLAR_DIRECTIONS) where {N}
    n = dim(cup)
    # extract limits and extrema of already plotted sets
    p = plotattributes[:plot_object]
    lims = _extract_limits(p, plotattributes)
    extr = _extract_extrema(p)
    cup_bounds = cup
    if isbounded(cup)
        B = box_approximation(cup)
        cup_bounds = B
        extr[:x] = (min(extr[:x][1], low(B, 1)), max(extr[:x][2], high(B, 1)))
        if n > 1
            extr[:y] = (min(extr[:y][1], low(B, 2)), max(extr[:y][2], high(B, 2)))
        end
    end
    # if there is already a plotted set and the limits are fixed,
    # automatically adjust the axis limits (e.g. after plotting an
    # unbounded set)
    _set_auto_limits_to_extrema!(lims, extr)
    cup_bounds = intersection(cup_bounds, _bounding_hyperrectangle(lims, n, eltype(cup_bounds)))
    _update_plot_limits!(lims, cup_bounds)
    xlims --> lims[:x]
    if n > 1
        ylims --> lims[:y]
    end

    if same_recipe
        label --> DEFAULT_LABEL
        grid --> DEFAULT_GRID
        if DEFAULT_ASPECT_RATIO != :none
            aspect_ratio --> DEFAULT_ASPECT_RATIO
        end
        seriesalpha --> DEFAULT_ALPHA
        seriescolor --> DEFAULT_COLOR
        seriestype --> :shape
        return _plot_list_same_recipe(array(cup), ε, Nφ)
    else
        for Xi in array(cup)
            if Xi isa Intersection
                @series Xi, ε, Nφ
            else
                @series Xi, ε
            end
        end
    end
end

@recipe function plot_polyzono(P::AbstractSparsePolynomialZonotope{N}, ε::Real=N(PLOT_PRECISION);  # COV_EXCL_LINE
                               nsdiv=10, partition=nothing) where {N}
    if iszero(ngens_dep(P))  # optimization
        return convert(Zonotope, P)
    else
        label --> DEFAULT_LABEL
        grid --> DEFAULT_GRID
        if DEFAULT_ASPECT_RATIO != :none
            aspect_ratio --> DEFAULT_ASPECT_RATIO
        end
        seriesalpha --> DEFAULT_ALPHA
        seriescolor --> DEFAULT_COLOR
        seriestype --> :shape
        Poa = overapproximate(P, UnionSetArray{Zonotope}; nsdiv, partition)
        return _plot_list_same_recipe(array(Poa), ε)
    end
end

#######################
# plot_recipe methods #
#######################

"""
    plot_recipe(X::LazySet, [ε])

Convert a compact convex set to a pair `(x, y)` of points for plotting.

### Input

- `X` -- compact convex set
- `ε` -- approximation-error bound

### Output

A pair `(x, y)` of points that can be plotted.

### Notes

We do not support four-dimensional or higher-dimensional sets at the moment.

### Algorithm

One-dimensional sets are converted to an `Interval`.

For two-dimensional sets, we first compute a polygonal overapproximation.
The second argument, `ε`, corresponds to the error in Hausdorff distance between
the overapproximating set and `X`.
On the other hand, if you only want to produce a fast box-overapproximation of
`X`, pass `ε=Inf`.

Finally, we use the plot recipe for the constructed set (interval or polygon).

For three-dimensional sets, we assume that `X` is a polytope.
"""
function plot_recipe(X::LazySet, ε)
    @assert dim(X) <= 3 "this implementation cannot plot $(dim(X))-dimensional sets"
    @assert isbounded(X) "this implementation can only plot bounded sets"
    @assert isconvex(X) "this implementation can only plot convex sets"

    if dim(X) == 1
        Y = convert(Interval, X)
    elseif dim(X) == 2
        Y = overapproximate(X, ε)
    else
        return _plot_recipe_3d_polytope(X)
    end
    return plot_recipe(Y, ε)
end

"""
    plot_recipe(P::AbstractPolyhedron{N}, [ε]=zero(N)) where {N}

Convert a (bounded) polyhedron to a pair `(x, y)` of points for plotting.

### Input

- `P` -- bounded polyhedron
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

A pair `(x, y)` of points that can be plotted, where `x` is the vector of
x-coordinates and `y` is the vector of y-coordinates.

### Algorithm

We first assert that `P` is bounded (i.e., that `P` is a polytope).

One-dimensional polytopes are converted to an `Interval`.
Three-dimensional or higher-dimensional polytopes are not supported.

For two-dimensional polytopes (i.e., polygons) we compute their set of vertices
using `vertices_list` and then plot the convex hull of these vertices.
"""
function plot_recipe(P::AbstractPolyhedron{N}, ε=zero(N)) where {N}
    @assert dim(P) <= 3 "cannot plot a $(dim(P))-dimensional $(typeof(P))"
    @assert isbounded(P) "cannot plot an unbounded $(typeof(P))"

    if dim(P) == 1
        Q = convert(Interval, P)
        if diameter(Q) < _ztol(N)  # flat interval
            Q = Singleton(center(Q))
        end
        return plot_recipe(Q, ε)
    elseif dim(P) == 2
        vlist = convex_hull(vertices_list(P))
        return _plot_recipe_2d_vlist(vlist, N)
    else
        return _plot_recipe_3d_polytope(P, N)
    end
end

function _plot_recipe_2d_vlist(vlist, N)
    m = length(vlist)
    if m == 0
        @warn "received a polyhedron with no vertices during plotting"
        return plot_recipe(EmptySet{N}(2), zero(N))
    end

    x = Vector{N}(undef, m)
    y = Vector{N}(undef, m)
    @inbounds for (i, vi) in enumerate(vlist)
        x[i] = vi[1]
        y[i] = vi[2]
    end

    if m > 2
        # add first vertex to "close" the polygon
        push!(x, x[1])
        push!(y, y[1])
    end
    return x, y
end

"""
    plot_recipe(S::AbstractSingleton{N}, [ε]=zero(N)) where {N}

Convert a singleton to a pair `(x, y)` of points for plotting.

### Input

- `S` -- singleton
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

A pair `(x, y)` of one point that can be plotted.
"""
function plot_recipe(S::AbstractSingleton{N}, ε=zero(N)) where {N}
    n = dim(S)
    if n == 1
        return [center(S)[1]], [zero(N)]
    elseif n == 2
        return [center(S)[1]], [center(S)[2]]
    elseif n == 3
        return [center(S)[1]], [center(S)[2]], [center(S)[3]]
    else
        throw(ArgumentError("cannot plot a $n-dimensional $(typeof(S))"))
    end
end

"""
    plot_recipe(cap::Intersection{N}, [ε]::N=-one(N),
                [Nφ]::Int=PLOT_POLAR_DIRECTIONS) where {N}

Convert an intersection of two sets to a pair `(x, y)` of points for plotting.

### Input

- `cap` -- intersection of two sets
- `ε`   -- (optional, default `0`) ignored, used for dispatch
- `Nφ`  -- (optional, default: `PLOT_POLAR_DIRECTIONS`) number of polar
           directions used in the template overapproximation

### Output

A pair `(x, y)` of points that can be plotted.
"""
function plot_recipe(cap::Intersection{N}, ε::N=zero(N),
                     Nφ::Int=PLOT_POLAR_DIRECTIONS) where {N}
    @assert dim(cap) <= 2 "cannot plot a $(dim(cap))-dimensional intersection"

    if isempty(cap)
        return plot_recipe(EmptySet{N}(dim(cap)), ε)
    elseif dim(cap) == 1
        @assert isconvex(cap) "cannot plot a one-dimensional $(typeof(cap))"
        return plot_recipe(convert(Interval, cap), ε)
    else
        # construct polygon approximation using polar directions
        P = overapproximate(cap, PolarDirections{N}(Nφ))
        return plot_recipe(P, ε)
    end
end

"""
    plot_recipe(∅::EmptySet{N}, [ε]=zero(N)) where {N}

Convert an empty set to a sequence of points for plotting.
In the special case of an empty set, the sequence is empty.

### Input

- `∅` -- empty set
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

An empty array.
"""
function plot_recipe(::EmptySet{N}, ε=zero(N)) where {N}
    return []
end

"""
    plot_recipe(X::Interval{N}, [ε]=zero(N)) where {N}

Convert an interval to a pair `(x, y)` of points for plotting.

### Input

- `X` -- interval
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

A pair `(x, y)` of two points that can be plotted.

### Notes

We consider the interval as a line segment with y coordinate equal to zero.
"""
function plot_recipe(X::Interval{N}, ε=zero(N)) where {N}
    bounds = extrema(X, 1)
    x = @inbounds [bounds[1], bounds[2]]
    return x, zeros(N, 2)
end

function plot_recipe(P::Polygon{N}, ε=zero(N)) where {N}
    return _plot_recipe_2d_vlist(P.vertices, N)
end

end  # module
