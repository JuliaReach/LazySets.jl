export ρ_upper_bound

"""
    ρ_upper_bound(d::AbstractVector{N},
                  X::LazySet) where {N<:Real}

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
                       X::LazySet) where {N<:Real}
    return ρ(d, X)
end

"""
    ρ_upper_bound(d::AbstractVector{N},
                  cap::Intersection{N, <:LazySet, S}) where {N<:Real, S<:AbstractPolytope{N}}

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
                       cap::Intersection{N, <:LazySet, S}) where {N<:Real, S<:AbstractPolytope{N}}
    return min(ρ_upper_bound(d, cap.X), ρ_upper_bound(d, cap.Y))
end

"""
    ρ_upper_bound(d::AbstractVector{N},
                  cap::Intersection{N, <:LazySet, S};
                  kwargs...) where {N<:Real, S<:AbstractPolytope{N}}

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
                       cap::Intersection{N, <:LazySet, S};
                       kwargs...) where {N<:Real, S<:AbstractPolytope{N}}

    X = cap.X    # compact set
    P = cap.Y    # polytope
    return min([ρ_upper_bound(d, X ∩ Hi, kwargs...) for Hi in constraints_list(P)])
end
