import Base.<=

export HPolygon

"""
    HPolygon{N<:Real} <: AbstractHPolygon{N}

Type that represents a convex polygon in constraint representation whose edges
are sorted in counter-clockwise fashion with respect to their normal directions.

### Fields

- `constraints`          -- list of linear constraints, sorted by the angle
- `sort_constraints`     -- (optional, default: `true`) flag for sorting the
                            constraints (sortedness is a running assumption of
                            this type)
- `validate_boundedness` -- (optional, default: `true`) flag for checking if the
                            constraints make the polygon bounded; (boundedness
                            is a running assumption of this type)

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

    # default constructor that applies sorting of the given constraints
    function HPolygon{N}(constraints::Vector{LinearConstraint{N}};
                         sort_constraints::Bool=true,
                         validate_boundedness::Bool=false) where {N<:Real}
        if sort_constraints
            sorted_constraints = Vector{LinearConstraint{N}}()
            sizehint!(sorted_constraints, length(constraints))
            for ci in constraints
                addconstraint!(sorted_constraints, ci)
            end
            P = new{N}(sorted_constraints)
        else
            P = new{N}(constraints)
        end
        @assert (!validate_boundedness ||
                 LazySets.validate_boundedness(P)) "the polygon is not bounded"
        return P
    end
end

# convenience constructor without type parameter
HPolygon(constraints::Vector{LinearConstraint{N}};
         sort_constraints::Bool=true,
         validate_boundedness::Bool=false) where {N<:Real} =
    HPolygon{N}(constraints;
                sort_constraints=sort_constraints,
                validate_boundedness=validate_boundedness)

# constructor for an HPolygon with no constraints
HPolygon{N}() where {N<:Real} = HPolygon{N}(Vector{LinearConstraint{N}}())

# constructor for an HPolygon with no constraints of type Float64
HPolygon() = HPolygon{Float64}()

# conversion constructor
function HPolygon(S::LazySet; validate_boundedness::Bool=false)
    P = convert(HPolygon, S)
    if validate_boundedness
        # trigger boundedness check in constructor
        HPolygon(P.constraints; validate_boundedness=true)
    end
    return P
end


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
        while k <= n && P.constraints[k].a <= d
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
