import Base.<=

export HPolygon

"""
    HPolygon{N<:Real} <: AbstractHPolygon{N}

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
struct HPolygon{N<:Real} <: AbstractHPolygon{N}
    constraints_list::Vector{LinearConstraint{N}}
end
# constructor for an HPolygon with no constraints
HPolygon{N}() where {N<:Real} = HPolygon{N}(Vector{N}(0))
# constructor for an HPolygon with no constraints of type Float64
HPolygon() = HPolygon{Float64}()


# --- LazySet interface functions ---


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
