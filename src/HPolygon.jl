import Base.<=

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
HPolygon{N}() where {N<:Real} = HPolygon{N}(Vector{LinearConstraint{N}}(0))

# constructor for an HPolygon with no constraints of type Float64
HPolygon() = HPolygon{Float64}()

# conversion constructor
HPolygon(S::LazySet) = convert(HPolygon, S)


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
    n = length(P.constraints)
    @assert n > 0 "the polygon has no constraints"
    k = 1
    while k <= n && P.constraints[k].a <= d
        k += 1
    end
    if k == 1 || k == n+1
        return element(intersection(Line(P.constraints[1]),
                                    Line(P.constraints[n])))
    else
        return element(intersection(Line(P.constraints[k]),
                                    Line(P.constraints[k-1])))
    end
end
