import Base.<=

export HPolygon

"""
    HPolygon{N<:Real} <: AbstractHPolygon{N}

Type that represents a convex polygon in constraint representation whose edges
are sorted in counter-clockwise fashion with respect to their normal directions.

### Fields

- `constraints`       -- list of linear constraints, sorted by the angle
- `sort_constraints`  -- (optional, default: `true`) flag for sorting the
                         constraints (sortedness is a running assumption of this
                         type)
- `check_boundedness` -- (optional, default: `false`) flag for checking if the
                         constraints make the polygon bounded; (boundedness is a
                         running assumption of this type)
- `prune`             -- (optional, default: `true`) flag for removing redundant
                         constraints

### Notes

The default constructor assumes that the given list of edges is sorted.
It *does not perform* any sorting.
Use `addconstraint!` to iteratively add the edges in a sorted way.

- `HPolygon(constraints::Vector{LinearConstraint{<:Real}})`
  -- default constructor
- `HPolygon()`
  -- constructor with no constraints
"""
struct HPolygon{N<:Real} <: AbstractHPolygon{N}
    constraints::Vector{LinearConstraint{N}}

    # default constructor that applies sorting of the given constraints and
    # (checks for and) removes redundant constraints
    function HPolygon{N}(constraints::Vector{LinearConstraint{N}};
                         sort_constraints::Bool=true,
                         check_boundedness::Bool=false,
                         prune::Bool=true) where {N<:Real}
        if sort_constraints
            sorted_constraints = Vector{LinearConstraint{N}}()
            sizehint!(sorted_constraints, length(constraints))
            for ci in constraints
                addconstraint!(sorted_constraints, ci; prune=prune)
            end
            P = new{N}(sorted_constraints)
        else
            P = new{N}(constraints)
        end
        @assert (!check_boundedness ||
                 isbounded(P, false)) "the polygon is not bounded"
        return P
    end
end

# convenience constructor without type parameter
HPolygon(constraints::Vector{LinearConstraint{N}};
         sort_constraints::Bool=true,
         check_boundedness::Bool=false,
         prune::Bool=true) where {N<:Real} =
    HPolygon{N}(constraints;
                sort_constraints=sort_constraints,
                check_boundedness=check_boundedness,
                prune=prune)

# constructor for an HPolygon with no constraints
HPolygon{N}() where {N<:Real} = HPolygon{N}(Vector{LinearConstraint{N}}())

# constructor for an HPolygon with no constraints of type Float64
HPolygon() = HPolygon{Float64}()

# constructor from a simple H-representation
HPolygon(A::AbstractMatrix{N},
         b::AbstractVector{N};
         sort_constraints::Bool=true,
         check_boundedness::Bool=false,
         prune::Bool=true) where {N<:Real} =
    HPolygon(constraints_list(A, b); sort_constraints=sort_constraints,
             check_boundedness=check_boundedness, prune=prune)

# constructor from a simple H-representation with type parameter
HPolygon{N}(A::AbstractMatrix{N},
            b::AbstractVector{N};
            sort_constraints::Bool=true,
            check_boundedness::Bool=false,
            prune::Bool=true) where {N<:Real} =
    HPolygon(A, b; sort_constraints=sort_constraints,
             check_boundedness=check_boundedness, prune=prune)


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

"""
    translate(v::AbstractVector{N}, P::HPolygon{N}; share::Bool=false
             ) where {N<:Real}

Translate (i.e., shift) a polygon in constraint representation by a given
vector.

### Input

- `P`     -- polygon in constraint representation
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated polygon in constraint representation.

### Notes

The normal vectors of the constraints (vector `a` in `a⋅x ≤ b`) are shared with
the original constraints if `share == true`.

### Algorithm

We translate every constraint.
"""
function translate(P::HPolygon{N}, v::AbstractVector{N}; share::Bool=false
                  ) where {N<:Real}
    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    constraints = [translate(c, v; share=share) for c in constraints_list(P)]
    return HPolygon(constraints;
                    sort_constraints=false, check_boundedness=false,
                    prune=false)
end
