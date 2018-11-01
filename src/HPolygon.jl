export HPolygon

"""
    HPolygon{N<:Real} <: AbstractHPolygon{N}

Type that represents a convex polygon in constraint representation whose edges
are sorted in counter-clockwise fashion with respect to their normal directions.

### Fields

- `constraints` -- list of linear constraints, sorted by the angle

### Notes

The default constructor assumes that the given list of edges is sorted.
It *does not perform* any sorting.
Use `addconstraint!` to iteratively add the edges in a sorted way.

- `HPolygon(constraints::Vector{LinearConstraint{<:Real}})`
  -- default constructor
- `HPolygon()`
  -- constructor with no constraints
- `HPolygon(S::LazySet)` -- constructor from another set
"""
struct HPolygon{N<:Real} <: AbstractHPolygon{N}
    constraints::Vector{LinearConstraint{N}}
end

# constructor for an HPolygon with no constraints
HPolygon{N}() where {N<:Real} = HPolygon{N}(Vector{LinearConstraint{N}}())

# constructor for an HPolygon with no constraints of type Float64
HPolygon() = HPolygon{Float64}()

# conversion constructor
HPolygon(S::LazySet) = convert(HPolygon, S)


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, P::HPolygon{N};
      [linear_search]::Bool=(length(P.constraints) < BINARY_SEARCH_THRESHOLD)
     ) where {N<:Real}

Return the support vector of a polygon in a given direction.

### Input

- `d`             -- direction
- `P`             -- polygon in constraint representation
- `linear_search` -- (optional, default: see below) flag for controlling whether
                     to perform a linear search or a binary search

### Output

The support vector in the given direction.
The result is always one of the vertices; in particular, if the direction has
norm zero, any vertex is returned.

### Algorithm

Comparison of directions is performed using polar angles; see the overload of
`<=` for two-dimensional vectors.

For polygons with `BINARY_SEARCH_THRESHOLD = 10` or more constraints we use a
binary search by default.
"""
function σ(d::AbstractVector{N}, P::HPolygon{N};
           linear_search::Bool=(length(P.constraints) < BINARY_SEARCH_THRESHOLD)
          ) where {N<:Real}
    n = length(P.constraints)
    @assert n > 0 "the polygon has no constraints"

    if linear_search
        # linear search
        k = 1
        quadrant_d = quadrant(d)
        while k <= n && le_quadrant(P.constraints[k].a, d, quadrant_d)
            k += 1
        end
    else
        # binary search
        k = binary_search_constraints(d, P.constraints, n, 1 + div(n, 2))
    end

    if k == 1 || k == n+1
        # corner cases: wrap-around in constraints list
        return element(intersection(Line(P.constraints[1]),
                                    Line(P.constraints[n])))
    else
        return element(intersection(Line(P.constraints[k]),
                                    Line(P.constraints[k-1])))
    end
end
