import Base.<=

# Polygon in constraint representation (optimized for large number of constraints)
export HPolygonOpt, addconstraint!, is_contained, tovrep, vertices_list

"""
    HPolygonOpt <: LazySet

Type that represents a convex polygon in constraint representation, whose edges
are sorted in counter-clockwise fashion with respect to their normal directions.
This is a refined version of `HPolygon`.

### Fields

- `P`   -- polygon
- `ind` -- an index in the list of constraints to begin the search to
           evaluate the support functions.

### Notes

This structure is optimized to evaluate the support function/vector with a large
sequence of directions, which are one to one close. The strategy is to have an
index that can be used to warm-start the search for optimal values in the support
vector computation.
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

- `P` -- optimized polygon
"""
function dim(P::HPolygonOpt)::Int64
    return 2
end

"""
    σ(d, P)

Return the support vector of the optimized polygon in a given direction.

### Input

- `d` -- direction
- `P` -- optimized polygon in constraint representation
"""
function σ(d::AbstractVector{Float64}, P::HPolygonOpt)::Vector{Float64}
    m = length(P.constraints_list)
    cl = P.constraints_list
    if (d <= cl[P.ind].a)
        k = P.ind-1
        while (k >= 1 && d <= cl[k].a)
            k -= 1
        end
        if (k == 0)
            P.ind = m
            return intersection(Line(cl[m]), Line(cl[1]))
        else
            P.ind = k
            return intersection(Line(cl[k]), Line(cl[k+1]))
        end
    else
        k = P.ind+1
        while (k <= m && cl[k].a <= d)
            k += 1
        end
        if (k == m+1)
            P.ind = m
            return intersection(Line(cl[m]), Line(cl[1]))
        else
            P.ind = k-1
            return intersection(Line(cl[k-1]), Line(cl[k]))
        end
    end
end

"""
    is_contained(x, P)

Return whether a given vector is contained in an optimized polygon in constraint
representation.

### Input

- `x` -- two-dimensional vector
- `P` -- optimized polygon in constraint representation

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

- `P` -- optimized polygon in constraint representation

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
