export ρ_upper_bound

"""
    ρ_upper_bound(d::AbstractVector{N},
                  X::LazySet{N};
                  kwargs...) where {N<:Real}

Return an upper bound of the support function of a given set.

### Input

- `d` -- direction
- `X` -- set

### Output

An upper bound of the support function of the given set.

### Algorithm

The default implementation of `ρ_upper_bound` is the exact `ρ(d, X)`.
"""
function ρ_upper_bound(d::AbstractVector{N},
                       X::LazySet{N};
                       kwargs...) where {N<:Real}
    return ρ(d, X)
end

"""
    ρ_upper_bound(d::AbstractVector{N},
                  X::LazySet{N};
                  kwargs...) where {N<:Real}

Return an upper bound of the support function of a given set.

### Input

- `d` -- direction
- `X` -- set

### Output

An upper bound of the support function of the given set.

### Algorithm

The default implementation of `ρ_upper_bound` for lazy operations is to pass
all keyword arguments to `ρ(d, X; kwargs...)` if that method exists.
"""
function ρ_upper_bound(d::AbstractVector{N},
                       X::Union{LinearMap{N}};
                       kwargs...) where {N<:Real}
    return ρ(d, X, upper_bound=true, kwargs...)
end

"""
    ρ_upper_bound(d::AbstractVector{N},
                  cap::Intersection{N}; kwargs...) where {N<:Real}

Return an upper bound of the support function of the intersection of two sets.

### Input

- `d`   -- direction
- `cap` -- intersection

### Output

An upper bound of the support function of the given intersection.

### Algorithm

The support function of an intersection of ``X`` and ``Y`` is upper bounded by
the minimum of the support functions of ``X`` and ``Y``.
"""
function ρ_upper_bound(d::AbstractVector{N},
                       cap::Intersection{N}; kwargs...) where {N<:Real}
    return min(ρ_upper_bound(d, cap.X; kwargs...), ρ_upper_bound(d, cap.Y; kwargs...))
end

"""
    ρ_upper_bound(d::AbstractVector{N},
                  cap::Intersection{N, <:LazySet{N}, <:AbstractPolytope{N}};
                  kwargs...) where {N<:Real}

Return an upper bound of the intersection between a compact set and a
polytope along a given direction.

### Input

- `d`      -- direction
- `cap`    -- intersection of a compact set and a polytope
- `kwargs` -- additional arguments that are passed to the support function algorithm

### Output

An upper bound of the support function of the given intersection.

### Algorithm

The idea is to solve the univariate optimization problem `ρ(di, X ∩ Hi)` for each
half-space in the set `P` and then take the minimum. This gives an overapproximation
of the exact support function.

This algorithm is inspired from [G. Frehse, R. Ray. Flowpipe-Guard Intersection
for Reachability Computations with Support
Functions](https://www.sciencedirect.com/science/article/pii/S1474667015371809).

### Notes

This method relies on having available the `constraints_list` of the polytope
`P`.

This method of overapproximation can return a non-empty set even if the original
intersection is empty.
"""
function ρ_upper_bound(d::AbstractVector{N},
                       cap::Intersection{N, <:LazySet{N}, <:AbstractPolytope{N}};
                       kwargs...) where {N<:Real}

    X = cap.X    # compact set
    P = cap.Y    # polytope
    return minimum([ρ(d, X ∩ Hi; kwargs...) for Hi in constraints_list(P)])
end

# disambiguation
function ρ_upper_bound(d::AbstractVector{N},
                       cap::Intersection{N, <:AbstractPolytope{N}, <:AbstractPolytope{N}};
                       kwargs...) where {N<:Real}
    X = cap.X    # compact set
    P = cap.Y    # polytope
    return minimum([ρ(d, X ∩ Hi; kwargs...) for Hi in constraints_list(P)])
end

# symmetric function
function ρ_upper_bound(d::AbstractVector{N},
                       cap::Intersection{N, <:AbstractPolytope{N}, <:LazySet{N}};
                       kwargs...) where {N<:Real}
    return ρ_upper_bound(d, cap.Y ∩ cap.X; kwargs...)
end

"""
    ρ_upper_bound(d::AbstractVector{N},
                  cap::Intersection{N, <:LazySet{N},
                                       <:Union{HalfSpace{N}, Hyperplane{N}, Line{N}}};
                  [algorithm]::String="line_search",
                  [check_intersection]::Bool=true,
                  [kwargs...]) where N<:Real

Return an upper bound of the support function of the intersection of a compact
set and a half-space or a hyperplane or a line in a given direction.

### Input

- `d`         -- direction
- `cap`       -- lazy intersection of a compact set and a half-space/hyperplane/
                 line
- `algorithm` -- (optional, default: `"line_search"`): the algorithm to
                 calculate the support function; valid options are:

    * `"line_search"` -- solve the associated univariate optimization problem
                         using a line search method (either Brent or the
                         Golden Section method)
    * `"projection"`  -- only valid for intersection with a hyperplane;
                         evaluates the support function by reducing the problem
                         to the 2D intersection of a rank 2 linear
                         transformation of the given compact set in the plane
                         generated by the given direction `d` and the
                         hyperplane's normal vector `n`

- `check_intersection` -- (optional, default: `true`) if `true`, check if the
                          intersection is empty before actually calculating the
                          support function

### Output

The scalar value of the support function of the set `cap` in the given
direction.

### Notes

It is assumed that the set `cap.X` is compact.

The `check_intersection` flag can be useful if it is known in advance that the
intersection is non-empty.

Any additional number of arguments to the algorithm backend can be passed as
keyword arguments.

### Algorithm

We fall back to the implementation of `ρ` with the same arguments.
"""
function ρ_upper_bound(d::AbstractVector{N},
                       cap::Intersection{N,
                                         <:LazySet{N},
                                         <:Union{HalfSpace{N}, Hyperplane{N}, Line{N}}};
                       algorithm::String="line_search",
                       check_intersection::Bool=true,
                       kwargs...
                      ) where N<:Real
    return ρ(d, cap; algorithm=algorithm, check_intersection=check_intersection,
             kwargs...)
end

# symmetric case
function ρ_upper_bound(d::AbstractVector{N},
                       cap::Intersection{N,
                                         <:Union{HalfSpace{N}, Hyperplane{N}, Line{N}},
                                         <:LazySet{N}};
                       algorithm::String="line_search",
                       check_intersection::Bool=true,
                       kwargs...
                      ) where N<:Real
    return ρ(d, cap; algorithm=algorithm, check_intersection=check_intersection,
             kwargs...)
end

# disambiguation
function ρ_upper_bound(d::AbstractVector{N},
                       cap::Intersection{N,
                                         <:Union{HalfSpace{N}, Hyperplane{N}, Line{N}},
                                         <:AbstractPolytope{N}};
                       algorithm::String="line_search",
                       check_intersection::Bool=true,
                       kwargs...
                      ) where N<:Real
    return ρ(d, cap; algorithm=algorithm, check_intersection=check_intersection,
             kwargs...)
end

# disambiguation symmetric case
function ρ_upper_bound(d::AbstractVector{N},
                       cap::Intersection{N,
                                         <:AbstractPolytope{N},
                                         <:Union{HalfSpace{N}, Hyperplane{N}, Line{N}}};
                       algorithm::String="line_search",
                       check_intersection::Bool=true,
                       kwargs...
                      ) where N<:Real
    return ρ(d, cap; algorithm=algorithm, check_intersection=check_intersection,
             kwargs...)
end
