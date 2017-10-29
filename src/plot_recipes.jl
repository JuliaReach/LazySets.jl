# ====================================
# Plot recipes for an abstract LazySet
# ====================================

"""
    plot_Polygon(P::Union{HPolygon, HPolygonOpt})

Plot a polygon given in constraint form.
 
### Input

- `X` -- a convex set

### Examples

```julia
julia> using LazySets, Plots
julia> B = BallInf(ones(2), 0.1)
julia> plot(2.0 * B)
```
"""
@recipe function plot_LazySet(X::T; color="blue", label="",
                              grid=true, alpha=0.5) where {T<:LazySet}

    seriestype := :shape

    P = Approximations.overapproximate(X)
    vlist = hcat(vertices_list(P)...).'
    (x, y) = vlist[:, 1], vlist[:, 2]
 
     x, y
end

"""
    plot_Polygon(P::Union{Vector{HPolygon}, Vector{HPolygonOpt}})

Plot an array of polygons given in constraint form.
 
### Input

- `P` -- an array of convex set

### Examples

```julia
julia> using LazySets, Plots
julia> B1 = BallInf(zeros(2), 0.4)
julia> B2 = BallInf(ones(2), 0.4)
julia> plot([B1, B2])
```
"""
@recipe function plot_LazySet(X::Vector{T}; color="blue", label="",
                              grid=true, alpha=0.5) where {T<:LazySet}

    seriestype := :shape

    for Xi in X
        Pi = Approximations.overapproximate(Xi)
        vlist = hcat(vertices_list(Pi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end
