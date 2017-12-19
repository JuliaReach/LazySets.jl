import Base.∈

export AbstractHPolygon,
       addconstraint!

"""
    AbstractHPolygon{N<:Real} <: AbstractPolygon{N}

Abstract type for polygons in H-representation (i.e., constraints).

### Notes

Every concrete `AbstractHPolygon` must have the following fields:
- `constraints_list::Vector{LinearConstraint{N}}` -- the constraints

```jldoctest
julia> subtypes(AbstractHPolygon)
2-element Array{Union{DataType, UnionAll},1}:
 LazySets.HPolygon   
 LazySets.HPolygonOpt
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
    vertices_list(P::AbstractHPolygon{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a polygon in constraint representation.

### Input

- `P` -- polygon in constraint representation

### Output

List of vertices.
"""
function vertices_list(P::AbstractHPolygon{N}
                      )::Vector{Vector{N}} where {N<:Real}
    n = length(P.constraints_list)
    points = Vector{Vector{N}}(n)
    if n == 0
        return points
    end
    @inbounds for i in 1:n-1
        points[i] = intersection(Line(P.constraints_list[i]),
                                 Line(P.constraints_list[i+1]))
    end
    points[n] = intersection(Line(P.constraints_list[n]),
                             Line(P.constraints_list[1]))
    return points
end


# --- LazySet interface functions ---


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

    for c in P.constraints_list
        if dot(c.a, x) > c.b
            return false
        end
    end
    return true
end


# --- common AbstractHPolygon functions ---


"""
    addconstraint!(P::AbstractHPolygon{N},
                   constraint::LinearConstraint{N})::Void where {N<:Real}

Add a linear constraint to a polygon in constraint representation, keeping the
constraints sorted by their normal directions.

### Input

- `P`          -- polygon in constraint representation
- `constraint` -- linear constraint to add

### Output

Nothing.
"""
function addconstraint!(P::AbstractHPolygon{N},
                        constraint::LinearConstraint{N})::Void where {N<:Real}
    i = length(P.constraints_list)
    while i > 0 && constraint.a <= P.constraints_list[i].a
        i -= 1
    end
    # here P.constraints_list[i] < constraint
    insert!(P.constraints_list, i+1, constraint)
    return nothing
end
