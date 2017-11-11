import Base.<=

# Polygon in constraint representation (optimized for large number of constraints)
export HPolygon, addconstraint!, is_contained, tovrep, vertices_list

"""
    HPolygonOpt <: LazySet

Type that represents a convex polygon (in H-representation).

### Fields

- `P`   -- polygon
- `ind` -- an index in the list of constraints to begin the search to
           evaluate the support functions.

### Notes

This structure is optimized to evaluate the support function/vector with a large
sequence of directions, which are one to one close.
"""
mutable struct HPolygonOpt <: LazySet
    constraints_list::Vector{LinearConstraint}
    ind::Int64

    HPolygonOpt(constraints_list) = new(constraints_list, 1)
    HPolygonOpt(constraints_list, ind) = new(constraints_list, ind)
end
HPolygonOpt(H::HPolygon) = HPolygonOpt(H.constraints_list)

"""
    dim(P)

Return the ambient dimension of the optimized polygon.

### Input

- `P` -- optimized polyhedron in H-representation
"""
function dim(P::HPolygonOpt)::Int64
    return 2
end

"""
    σ(d, P)

Return the support vector of the optimized polygon in a given direction.

### Input

- `d` -- direction
- `P` -- polyhedron in H-representation
"""
function σ(d::AbstractVector{Float64}, p::HPolygonOpt)::Vector{Float64}
    n = length(p.constraints_list)
    cl = p.constraints_list
    if (d <= cl[p.ind].a)
        k = p.ind-1
        while (k >= 1 && d <= cl[k].a)
            k -= 1
        end
        if (k == 0)
            p.ind = n
            return intersection(Line(cl[n]), Line(cl[1]))
        else
            p.ind = k
            return intersection(Line(cl[k]), Line(cl[k+1]))
        end
    else
        k = p.ind+1
        while (k <= n && cl[k].a <= d)
            k += 1
        end
        if (k == n+1)
            p.ind = n
            return intersection(Line(cl[n]), Line(cl[1]))
        else
            p.ind = k-1
            return intersection(Line(cl[k-1]), Line(cl[k]))
        end
    end
end

"""
    is_contained(x, P)

Return whether a given vector is contained in a polygon in constraint representation.

### Input

- `x` -- vector
- `P` -- polygon

### Output

Return `true` iff x ∈ P.
"""
function is_contained(x::AbstractVector{Float64}, P::HPolygonOpt)::Bool
    return is_contained(x, HPolygon(P.constraints_list))
end

"""
    tovrep(P)

Build a vertex representation of the given optimized polygon.

### Input

- `P` -- an optimized polygon in constraint representation

### Output

The same polygon in vertex representation, `VPolygon`.

### Note

The linear constraints of the input `HPolygonOpt` are assumed to be sorted by their
normal directions in counter-clockwise fashion.
"""
function tovrep(P::HPolygonOpt)::VPolygon
    return tovrep(HPolygon(P.constraints_list))
end

"""
    vertices_list(P)

Return the list of vertices of a convex polygon in constraint representation.

### Input

- `P` -- an optimized polygon in constraint representation

### Output

List of vertices as an array of vertex pairs, `Vector{Vector{Float64}}`.
"""
function vertices_list(P::HPolygonOpt)::Vector{Vector{Float64}}
    return vertices_list(HPolygon(P.constraints_list))
end
