import Base.<=

# Polygon in constraint representation
export HPolygon, addconstraint!, is_contained, tovrep, vertices_list

"""
    HPolygon <: LazySet

Type that represents a convex polygon in constraint representation, whose edges
are sorted in counter-clockwise fashion with respect to their normal directions.

### Fields

- `constraints_list` --  an array of linear constraints

### Note

The `HPolygon` constructor *does not perform* sorting of the given list of edges.
Use `addconstraint!` to iteratively add and sort the edges.
"""
struct HPolygon <: LazySet
    constraints_list::Vector{LinearConstraint}
end
HPolygon() = HPolygon([])

"""
    addconstraint!(P, constraint)

Add a linear constraint to a polygon in constraint representation keeping the
constraints sorted by their normal directions.

### Input

- `P`          -- a polygon
- `constraint` -- the linear constraint to add, see `LinearConstraint`
"""
function addconstraint!(P::HPolygon, constraint::LinearConstraint)
    i = length(P.constraints_list)
    while i > 0 && constraint.a <= P.constraints_list[i].a
        i -= 1
    end
    # here P.constraints_list[i] < constraint
    insert!(P.constraints_list, i+1, constraint)
end

"""
    dim(P)

Return the ambient dimension of the polygon.

### Input

- `P` -- polygon in constraint representation
"""
function dim(P::HPolygon)
    2
end

"""
    σ(d, P)

Return the support vector of a polygon in a given direction.

### Input

- `d` -- direction
- `P` -- polygon in constraint representation

### Algorithm

Comparison of directions is performed using polar angles, see the overload of
`<=` for two-dimensional vectors.
"""
function σ(d::AbstractVector{Float64}, P::HPolygon)::Vector{Float64}
    n = length(P.constraints_list)
    if n == 0
        error("this polygon is empty")
    end
    i = 1
    while i <= n && P.constraints_list[i].a <= d
        i += 1
    end
    if i == 1 || i == n+1
        intersection(Line(P.constraints_list[1]), Line(P.constraints_list[n]))
    else
        intersection(Line(P.constraints_list[i]), Line(P.constraints_list[i-1]))
    end
end

"""
    is_contained(x, P)

Return whether a given vector is contained in the polygon.

### Input

- `x` -- two-dimensional vector
- `P` -- polygon in constraint representation

### Output

Return `true` iff x ∈ P.
"""
function is_contained(x::AbstractVector{Float64}, P::HPolygon)::Bool
    res = true
    for c in P.constraints_list
        res = res && dot(c.a, x) <= c.b
    end
    return res
end

"""
    tovrep(P)

Build a vertex representation of the given polygon.

### Input

- `P` -- polygon in constraint representation

### Output

The same polygon but in vertex representation, `VPolygon`.

### Note

The linear constraints of the input `HPolygon` are assumed to be sorted by their
normal directions in counter-clockwise fashion.
"""
function tovrep(P::HPolygon)::VPolygon
    return VPolygon(vertices_list(P))
end

"""
    vertices_list(P)

Return the list of vertices of a convex polygon in constraint representation.

### Input

- `P` -- polygon in constraint representation

### Output

List of vertices as an array of vertex pairs, `Vector{Vector{Float64}}`.
"""
function vertices_list(P::HPolygon)::Vector{Vector{Float64}}
    m = length(P.constraints_list)
    points = Vector{Vector{Float64}}(m)
    @inbounds for i in 1:m-1
        points[i] = intersection(Line(P.constraints_list[i]), Line(P.constraints_list[i+1]))
    end
    points[m] = intersection(Line(P.constraints_list[m]), Line(P.constraints_list[1]))
    return points
end
