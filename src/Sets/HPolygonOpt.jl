import Base.<=

export HPolygonOpt

"""
    HPolygonOpt{N<:Real, VN<:AbstractVector{N}} <: AbstractHPolygon{N}

Type that represents a convex polygon in constraint representation whose edges
are sorted in counter-clockwise fashion with respect to their normal directions.
This is a refined version of [`HPolygon`](@ref).

### Fields

- `constraints`       -- list of linear constraints, sorted by the normal
                         direction in counter-clockwise fashion
- `ind`               -- index in the list of constraints to begin the search
                         to evaluate the support function
- `sort_constraints`  -- (optional, default: `true`) flag for sorting the
                         constraints (sortedness is a running assumption of this
                         type)
- `check_boundedness` -- (optional, default: `false`) flag for checking if the
                         constraints make the polygon bounded; (boundedness is a
                         running assumption of this type)
- `prune`             -- (optional, default: `true`) flag for removing redundant
                         constraints

### Notes

This structure is optimized to evaluate the support function/vector with a large
sequence of directions that are close to each other. The strategy is to have an
index that can be used to warm-start the search for optimal values in the
support vector computation.

The option `sort_constraints` can be used to deactivate automatic sorting of
constraints in counter-clockwise fashion, which is an invariant of this type.
Alternatively, one can construct an `HPolygonOpt` with empty constraints list,
which can then be filled iteratively using `addconstraint!`.

Similarly, the option `prune` can be used to deactivate automatic pruning of
redundant constraints.

Another type assumption is that the polygon is bounded.
The option `check_boundedness` can be used to assert this.
This option is deactivated by default because we explicitly want to allow the
iterative addition of the constraints, and hence one has to initially construct
an empty list of constraints (which represents an unbounded set).
The user has to make sure that the `HPolygonOpt` is not used before the
constraints actually describe a bounded set.
The function `isbounded` can be used to manually assert boundedness.
"""
mutable struct HPolygonOpt{N<:Real, VN<:AbstractVector{N}} <: AbstractHPolygon{N}
    constraints::Vector{LinearConstraint{N, VN}}
    ind::Int

    # default constructor that applies sorting of the given constraints
    function HPolygonOpt(constraints::Vector{LinearConstraint{N, VN}},
                         ind::Int=1;
                         sort_constraints::Bool=true,
                         check_boundedness::Bool=false,
                         prune::Bool=true) where {N<:Real, VN<:AbstractVector{N}}
        if sort_constraints
            sorted_constraints = Vector{LinearConstraint{N, VN}}()
            sizehint!(sorted_constraints, length(constraints))
            for ci in constraints
                addconstraint!(sorted_constraints, ci; prune=prune)
            end
            P = new{N, VN}(sorted_constraints, ind)
        else
            P = new{N, VN}(constraints, ind)
        end
        @assert (!check_boundedness ||
                 isbounded(P, false)) "the polygon is not bounded"
        return P
    end
end

isoperationtype(::Type{<:HPolygonOpt}) = false
isconvextype(::Type{<:HPolygonOpt}) = true

# constructor for an HPolygonOpt with no constraints and given numeric type
function HPolygonOpt{N}() where {N<:Real}
    HPolygonOpt(Vector{LinearConstraint{N, Vector{N}}}())
end

# constructor for an HPolygonOpt without explicit numeric type, defaults to Float64
function HPolygonOpt() where {N<:Real}
    HPolygonOpt{Float64}()
end

# constructor from a simple H-representation
HPolygonOpt(A::AbstractMatrix{N},
            b::AbstractVector{N};
            sort_constraints::Bool=true,
            check_boundedness::Bool=false,
            prune::Bool=true) where {N<:Real} =
    HPolygonOpt(constraints_list(A, b); sort_constraints=sort_constraints,
                check_boundedness=check_boundedness, prune=prune)


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, P::HPolygonOpt{N};
      [linear_search]::Bool=(length(P.constraints) < BINARY_SEARCH_THRESHOLD)
     ) where {N<:Real}

Return the support vector of an optimized polygon in a given direction.

### Input

- `d`             -- direction
- `P`             -- optimized polygon in constraint representation
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
function σ(d::AbstractVector{N}, P::HPolygonOpt{N};
           linear_search::Bool=(length(P.constraints) < BINARY_SEARCH_THRESHOLD)
          ) where {N<:Real}
    n = length(P.constraints)
    @assert n > 0 "the polygon has no constraints"
    if linear_search
        # linear search
        if (d <= P.constraints[P.ind].a)
            # search backward
            k = P.ind-1
            while (k >= 1 && d <= P.constraints[k].a)
                k -= 1
            end
            if (k == 0)
                P.ind = n
                # corner case: wrap-around in constraints list
                return element(intersection(Line(P.constraints[n]),
                                            Line(P.constraints[1])))
            else
                P.ind = k
            end
        else
            # search forward
            k = P.ind+1
            while (k <= n && P.constraints[k].a <= d)
                k += 1
            end
            if (k == n+1)
                P.ind = n
                # corner case: wrap-around in constraints list
                return element(intersection(Line(P.constraints[n]),
                                            Line(P.constraints[1])))
            else
                P.ind = k-1
            end
        end
        return element(intersection(Line(P.constraints[P.ind]),
                                    Line(P.constraints[P.ind + 1])))
    else
        # binary search
        k = binary_search_constraints(d, P.constraints, n, P.ind)
        if k == 1 || k == n+1
            P.ind = 1
            # corner cases: wrap-around in constraints list
            return element(intersection(Line(P.constraints[n]),
                                        Line(P.constraints[1])))
        else
            P.ind = k
            return element(intersection(Line(P.constraints[k-1]),
                                        Line(P.constraints[k])))
        end
    end
end

"""
    translate(P::HPolygonOpt{N}, v::AbstractVector{N}; share::Bool=false
             ) where {N<:Real}

Translate (i.e., shift) an optimized polygon in constraint representation by a
given vector.

### Input

- `P`     -- optimized polygon in constraint representation
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated optimized polygon in constraint representation.

### Notes

The normal vectors of the constraints (vector `a` in `a⋅x ≤ b`) are shared with
the original constraints if `share == true`.

### Algorithm

We translate every constraint.
"""
function translate(P::HPolygonOpt{N}, v::AbstractVector{N}; share::Bool=false
                  ) where {N<:Real}
    @assert length(v) == dim(P) "cannot translate a $(dim(P))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    constraints = [translate(c, v; share=share) for c in constraints_list(P)]
    return HPolygonOpt(constraints, P.ind;
                       sort_constraints=false, check_boundedness=false,
                       prune=false)
end
