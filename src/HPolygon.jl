import Base.<=

export HPolygon,
       addconstraint!,
       is_contained,
       tovrep,
       vertices_list

"""
    HPolygon{N<:Real} <: LazySet

Type that represents a convex polygon in constraint representation whose edges
are sorted in counter-clockwise fashion with respect to their normal directions.

### Fields

- `constraints_list` -- list of linear constraints, sorted by the angle

### Notes

The default constructor assumes that the given list of edges is sorted.
It *does not perform* any sorting.
Use `addconstraint!` to iteratively add the edges in a sorted way.

- `HPolygon(constraints_list::Vector{LinearConstraint{<:Real}})`
  -- default constructor
- `HPolygon()`
  -- constructor with no constraints
"""
struct HPolygon{N<:Real} <: LazySet
    constraints_list::Vector{LinearConstraint{N}}
end
# constructor for an HPolygon with no constraints
HPolygon{N}() where {N<:Real} = HPolygon{N}(Vector{N}(0))
# constructor for an HPolygon with no constraints of type Float64
HPolygon() = HPolygon{Float64}()

"""
    addconstraint!(P::HPolygon{N}, constraint::LinearConstraint{N}) where {N<:Real}

Add a linear constraint to a polygon in constraint representation, keeping the
constraints sorted by their normal directions.

### Input

- `P`          -- polygon
- `constraint` -- linear constraint to add
"""
function addconstraint!(P::HPolygon{N},
                        constraint::LinearConstraint{N}) where {N<:Real}
    i = length(P.constraints_list)
    while i > 0 && constraint.a <= P.constraints_list[i].a
        i -= 1
    end
    # here P.constraints_list[i] < constraint
    insert!(P.constraints_list, i+1, constraint)
end

"""
    dim(P::HPolygon)

Return the dimension of a polygon.

### Input

- `P` -- polygon in constraint representation

### Output

The ambient dimension of the polygon.
"""
function dim(P::HPolygon)
    return 2
end

"""
    σ(d::AbstractVector{<:Real}, P::HPolygon{N})::Vector{N} where {N<:Real}

Return the support vector of a polygon in a given direction.

### Input

- `d` -- direction
- `P` -- polygon in constraint representation

### Output

The support vector in the given direction.
The result is always one of the vertices; in particular, if the direction has
norm zero, any vertex is returned.

### Algorithm

Comparison of directions is performed using polar angles; see the overload of
`<=` for two-dimensional vectors.
"""
function σ(d::AbstractVector{<:Real}, P::HPolygon{N})::Vector{N} where {N<:Real}
    n = length(P.constraints_list)
    if n == 0
        error("this polygon is empty")
    end
    k = 1
    while k <= n && P.constraints_list[k].a <= d
        k += 1
    end
    if k == 1 || k == n+1
        return intersection(Line(P.constraints_list[1]),
                            Line(P.constraints_list[n]))
    else
        return intersection(Line(P.constraints_list[k]),
                            Line(P.constraints_list[k-1]))
    end
end

"""
    is_contained(x::AbstractVector{<:Real}, P::HPolygon)::Bool

Return whether a given vector is contained in a polygon.

### Input

- `x` -- two-dimensional vector
- `P` -- polygon in constraint representation

### Output

Return `true` iff ``x ∈ P``.
"""
function is_contained(x::AbstractVector{<:Real}, P::HPolygon)::Bool
    for c in P.constraints_list
        if !(dot(c.a, x) <= c.b)
            return false
        end
    end
    return true
end

"""
    tovrep(P::HPolygon)::VPolygon

Build a vertex representation of the given polygon.

### Input

- `P` -- polygon in constraint representation

### Output

The same polygon but in vertex representation, a `VPolygon`.
"""
function tovrep(P::HPolygon)::VPolygon
    return VPolygon(vertices_list(P))
end

"""
    vertices_list(P::HPolygon{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a polygon in constraint representation.

### Input

- `P` -- polygon in constraint representation

### Output

List of vertices.
"""
function vertices_list(P::HPolygon{N})::Vector{Vector{N}} where {N<:Real}
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
