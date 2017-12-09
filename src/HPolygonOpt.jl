import Base: <=, ∈

export HPolygonOpt,
       addconstraint!,
       ∈,
       tovrep,
       vertices_list

"""
    HPolygonOpt{N<:Real} <: LazySet

Type that represents a convex polygon in constraint representation whose edges
are sorted in counter-clockwise fashion with respect to their normal directions.
This is a refined version of `HPolygon`.

### Fields

- `constraints_list` -- list of linear constraints
- `ind` -- index in the list of constraints to begin the search to evaluate the
           support function

### Notes

This structure is optimized to evaluate the support function/vector with a large
sequence of directions that are close to each other. The strategy is to have an
index that can be used to warm-start the search for optimal values in the
support vector computation.

The default constructor assumes that the given list of edges is sorted.
It *does not perform* any sorting.
Use `addconstraint!` to iteratively add the edges in a sorted way.

- `HPolygonOpt(constraints_list::Vector{LinearConstraint{<:Real}}, ind::Int)`
  -- default constructor
- `HPolygonOpt(constraints_list::Vector{LinearConstraint{<:Real}})`
  -- constructor without index
- `HPolygonOpt(H::HPolygon{<:Real})`
  -- constructor from an HPolygon
"""
mutable struct HPolygonOpt{N<:Real} <: LazySet
    constraints_list::Vector{LinearConstraint{N}}
    ind::Int

    # default constructor
    HPolygonOpt{N}(constraints_list::Vector{LinearConstraint{N}},
                   ind::Int) where {N<:Real} =
        new{N}(constraints_list, ind)
end
# type-less convenience constructor
HPolygonOpt(constraints_list::Vector{LinearConstraint{N}},
            ind::Int) where {N<:Real} =
    HPolygonOpt{N}(constraints_list, ind)

# type-less convenience constructor without index
HPolygonOpt(constraints_list::Vector{LinearConstraint{N}}) where {N<:Real} =
    HPolygonOpt{N}(constraints_list, 1)

# constructor from an HPolygon
HPolygonOpt(H::HPolygon{N}) where {N<:Real} =
    HPolygonOpt{N}(H.constraints_list, 1)

"""
    addconstraint!(P::HPolygonOpt{N}, constraint::LinearConstraint{N})::Void where {N<:Real}

Add a linear constraint to an optimized polygon in constraint representation,
keeping the constraints sorted by their normal directions.

### Input

- `P`          -- optimized polygon
- `constraint` -- linear constraint to add

### Output

Nothing.
"""
function addconstraint!(P::HPolygonOpt{N},
                        constraint::LinearConstraint{N})::Void where {N<:Real}
    i = length(P.constraints_list)
    while i > 0 && constraint.a <= P.constraints_list[i].a
        i -= 1
    end
    insert!(P.constraints_list, i+1, constraint)
    return nothing
end

"""
    dim(P::HPolygonOpt)::Int

Return the dimension of an optimized polygon.

### Input

- `P` -- optimized polygon in constraint representation

### Output

The ambient dimension of the optimized polygon.
"""
function dim(P::HPolygonOpt)::Int
    return 2
end

"""
    σ(d::AbstractVector{<:Real}, P::HPolygonOpt{N})::Vector{N} where {N<:Real}

Return the support vector of an optimized polygon in a given direction.

### Input

- `d` -- direction
- `P` -- optimized polygon in constraint representation

### Output

The support vector in the given direction.
The result is always one of the vertices; in particular, if the direction has
norm zero, any vertex is returned.

### Algorithm

Comparison of directions is performed using polar angles; see the overload of
`<=` for two-dimensional vectors.
"""
function σ(d::AbstractVector{<:Real},
           P::HPolygonOpt{N})::Vector{N} where {N<:Real}
    n = length(P.constraints_list)
    if n == 0
        error("this polygon is empty")
    end
    if (d <= P.constraints_list[P.ind].a)
        k = P.ind-1
        while (k >= 1 && d <= P.constraints_list[k].a)
            k -= 1
        end
        if (k == 0)
            P.ind = n
            return intersection(Line(P.constraints_list[n]),
                                Line(P.constraints_list[1]))
        else
            P.ind = k
            return intersection(Line(P.constraints_list[k]),
                                Line(P.constraints_list[k+1]))
        end
    else
        k = P.ind+1
        while (k <= n && P.constraints_list[k].a <= d)
            k += 1
        end
        if (k == n+1)
            P.ind = n
            return intersection(Line(P.constraints_list[n]),
                                Line(P.constraints_list[1]))
        else
            P.ind = k-1
            return intersection(Line(P.constraints_list[k-1]),
                                Line(P.constraints_list[k]))
        end
    end
end

"""
    ∈(x::AbstractVector{<:Real}, P::HPolygonOpt)::Bool

Return whether a given vector is contained in an optimized polygon.

### Input

- `x` -- two-dimensional vector
- `P` -- optimized polygon in constraint representation

### Output

Return `true` iff ``x ∈ P``.
"""
function ∈(x::AbstractVector{<:Real}, P::HPolygonOpt)::Bool
    return ∈(x, HPolygon(P.constraints_list))
end

"""
    tovrep(P::HPolygonOpt)::VPolygon

Build a vertex representation of the given optimized polygon.

### Input

- `P` -- optimized polygon in constraint representation

### Output

The same polygon but in vertex representation, a `VPolygon`.
"""
function tovrep(P::HPolygonOpt)::VPolygon
    return tovrep(HPolygon(P.constraints_list))
end

"""
    vertices_list(P::HPolygonOpt{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of an optimized polygon in constraint
representation.

### Input

- `P` -- an optimized polygon in constraint representation

### Output

List of vertices.
"""
function vertices_list(P::HPolygonOpt{N})::Vector{Vector{N}} where {N<:Real}
    return vertices_list(HPolygon(P.constraints_list))
end
