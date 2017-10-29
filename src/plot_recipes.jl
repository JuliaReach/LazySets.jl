# ====================================
# Plot recipes for an abstract LazySet
# ====================================

"""
    plot_Polygon(X::T; ...) where {T<:LazySet}

Plot a lazy set in two-dimensions using an axis-aligned approximation.
 
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
    plot_Polygon(X::Vector{T}) where {T<:LazySet}

Plot an array of lazy sets in two-dimensions using an axis-aligned approximation.
 
### Input

- `X` -- a convex set

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

"""
    plot_Polygon(X::T, ε::Float64; ...) where {T<:LazySet}

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
@recipe function plot_LazySet(X::T, ε::Float64; color="blue", label="",
                              grid=true, alpha=0.5) where {T<:LazySet}

    seriestype := :shape

    P = Approximations.overapproximate(X, ε)
    vlist = hcat(vertices_list(P)...).'
    (x, y) = vlist[:, 1], vlist[:, 2]
 
     x, y
end

"""
    plot_Polygon(X::Vector{T}, ε::Float64; ...) where {T<:LazySet}

Plot an array of lazy sets in two-dimensions using iterative refinement.
 
### Input

- `X` -- a convex set
- `ε` -- approximation error bound

### Examples

```julia
julia> using LazySets, Plots
julia> B1 = BallInf(zeros(2), 0.4)
julia> B2 = Ball2(ones(2), 0.4)
julia> plot([B1, B2], 1e-4)
```
"""
@recipe function plot_LazySet(X::Vector{T}, ε::Float64; color="blue", label="",
                              grid=true, alpha=0.5) where {T<:LazySet}

    seriestype := :shape

    for Xi in X
        Pi = Approximations.overapproximate(Xi, ε)
        vlist = hcat(vertices_list(Pi)...).'
        @series (x, y) = vlist[:, 1], vlist[:, 2]
    end
end
