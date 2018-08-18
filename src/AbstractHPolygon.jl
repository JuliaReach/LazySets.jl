import Base.∈

export AbstractHPolygon,
       an_element,
       addconstraint!,
       vertices_list,
       constraints_list

# This constant marks the threshold for the number of constraints of a polygon
# above which we use a binary search to find the relevant constraint in a
# support vector query.
#
# NOTE: The value must be strictly greater than 2.
const BINARY_SEARCH_THRESHOLD = 10

"""
    AbstractHPolygon{N<:Real} <: AbstractPolygon{N}

Abstract type for polygons in H-representation (i.e., constraints).

### Notes

Every concrete `AbstractHPolygon` must have the following fields:
- `constraints::Vector{LinearConstraint{N}}` -- the constraints

New subtypes should be added to the `convert` method in order to be convertible.

```jldoctest
julia> subtypes(AbstractHPolygon)
2-element Array{Any,1}:
 HPolygon
 HPolygonOpt
```
"""
abstract type AbstractHPolygon{N<:Real} <: AbstractPolygon{N} end


# --- AbstractPolygon interface functions ---


"""
    tovrep(P::AbstractHPolygon{N})::VPolygon{N} where {N<:Real}

Build a vertex representation of the given polygon.

### Input

- `P` -- polygon in constraint representation

### Output

The same polygon but in vertex representation, a `VPolygon`.
"""
function tovrep(P::AbstractHPolygon{N})::VPolygon{N} where {N<:Real}
    return VPolygon(vertices_list(P))
end


"""
    tohrep(P::AbstractHPolygon{N})::AbstractHPolygon{N} where {N<:Real}

Build a contraint representation of the given polygon.

### Input

- `P` -- polygon in constraint representation

### Output

The identity, i.e., the same polygon instance.
"""
function tohrep(P::AbstractHPolygon{N})::AbstractHPolygon{N} where {N<:Real}
    return P
end


# --- AbstractPolytope interface functions ---


"""
    vertices_list(P::AbstractHPolygon{N},
                  apply_convex_hull::Bool=false
                 )::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a polygon in constraint representation.

### Input

- `P`                 -- polygon in constraint representation
- `apply_convex_hull` -- (optional, default: `false`) to post process or not the
                         intersection of constraints with a convex hull

### Output

List of vertices.
"""
function vertices_list(P::AbstractHPolygon{N},
                       apply_convex_hull::Bool=false
                      )::Vector{Vector{N}} where {N<:Real}
    n = length(P.constraints)
    points = Vector{Vector{N}}(undef, n)
    if n == 0
        return points
    end
    @inbounds for i in 1:n-1
        points[i] = element(intersection(Line(P.constraints[i]),
                                         Line(P.constraints[i+1])))
    end
    points[n] = element(intersection(Line(P.constraints[n]),
                                     Line(P.constraints[1])))
    return apply_convex_hull ? convex_hull(points) : points
end

"""
    constraints_list(P::AbstractHPolygon{N})::Vector{LinearConstraint{N}} where {N<:Real}

Return the list of constraints defining a polygon in H-representation.

### Input

- `P` -- polygon in H-representation

### Output

The list of constraints of the polygon.
"""
function constraints_list(P::AbstractHPolygon{N})::Vector{LinearConstraint{N}} where {N<:Real}
    return P.constraints
end

# --- LazySet interface functions ---


"""
    an_element(P::AbstractHPolygon{N})::Vector{N} where {N<:Real}

Return some element of a polygon in constraint representation.

### Input

- `P` -- polygon in constraint representation

### Output

A vertex of the polygon in constraint representation (the first one in the order
of the constraints).
"""
function an_element(P::AbstractHPolygon{N})::Vector{N} where {N<:Real}
    @assert length(P.constraints) >= 2 "polygon has less than two constraints"
    return element(intersection(Line(P.constraints[1]),
                                Line(P.constraints[2])))
end

"""
    ∈(x::AbstractVector{N}, P::AbstractHPolygon{N})::Bool where {N<:Real}

Check whether a given 2D point is contained in a polygon in constraint
representation.

### Input

- `x` -- two-dimensional point/vector
- `P` -- polygon in constraint representation

### Output

`true` iff ``x ∈ P``.

### Algorithm

This implementation checks if the point lies on the outside of each edge.
"""
function ∈(x::AbstractVector{N}, P::AbstractHPolygon{N})::Bool where {N<:Real}
    @assert length(x) == 2

    for c in P.constraints
        if dot(c.a, x) > c.b
            return false
        end
    end
    return true
end


# --- common AbstractHPolygon functions ---


"""
    addconstraint!(P::AbstractHPolygon{N},
                   constraint::LinearConstraint{N};
                   [linear_search]::Bool=(
                    length(P.constraints) < BINARY_SEARCH_THRESHOLD)
                  )::Nothing where {N<:Real}

Add a linear constraint to a polygon in constraint representation, keeping the
constraints sorted by their normal directions.

### Input

- `P`          -- polygon in constraint representation
- `constraint` -- linear constraint to add

### Output

Nothing.
"""
function addconstraint!(P::AbstractHPolygon{N},
                        constraint::LinearConstraint{N};
                        linear_search::Bool=(
                         length(P.constraints) < BINARY_SEARCH_THRESHOLD)
                       )::Nothing where {N<:Real}
    k = length(P.constraints)
    if k > 0
        d = constraint.a
        if d <= P.constraints[1].a
            k = 0
        else
            if linear_search
                # linear search
                while d <= P.constraints[k].a
                    k -= 1
                end
            else
                # binary search
                k = binary_search_constraints(
                    d, P.constraints, k, 1 + div(k, 2), choose_lower=true)
            end
        end
    end
    # here P.constraints[k] < constraint
    insert!(P.constraints, k+1, constraint)
    return nothing
end

"""
    binary_search_constraints(d::AbstractVector{N},
                              constraints::Vector{LinearConstraint{N}},
                              n::Int,
                              k::Int;
                              [choose_lower]::Bool=false
                             )::Int where {N<:Real}

Performs a binary search in the constraints.

### Input

- `d`            -- direction
- `constraints`  -- constraints
- `n`            -- number of constraints
- `k`            -- start index
- `choose_lower` -- (optional, default: `false`) flag for choosing the lower
                    index (see the 'Output' section)

### Output

In the default setting, the result is the smallest index `k` such that
`d <= constraints[k]`, or `n+1` if no such `k` exists.
If the `choose_lower` flag is set, the result is the largest index `k` such
that `constraints[k] < d`, which is equivalent to being `k-1` in the normal
setting.
"""
function binary_search_constraints(d::AbstractVector{N},
                                   constraints::Vector{LinearConstraint{N}},
                                   n::Int,
                                   k::Int;
                                   choose_lower::Bool=false
                                  )::Int where {N}
    lower = 1
    upper = n+1
    while lower + 1 < upper
        if constraints[k].a <= d
            lower = k
        else
            upper = k
        end
        k = lower + div(upper - lower, 2)
    end
    if choose_lower
        return lower
    else
        if lower == 1 && !(constraints[1].a <= d)
            # special case for index 1
            return 1
        end
        return upper
    end
end
