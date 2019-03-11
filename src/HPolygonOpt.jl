import Base.<=

export HPolygonOpt

"""
    HPolygonOpt{N<:Real} <: AbstractHPolygon{N}

Type that represents a convex polygon in constraint representation whose edges
are sorted in counter-clockwise fashion with respect to their normal directions.
This is a refined version of `HPolygon`.

### Fields

- `constraints`       -- list of linear constraints
- `ind`               -- index in the list of constraints to begin the search
                         to evaluate the support function
- `sort_constraints`  -- (optional, default: `true`) flag for sorting the
                         constraints (sortedness is a running assumption of this
                         type)
- `check_boundedness` -- (optional, default: `false`) flag for checking if the
                         constraints make the polygon bounded; (boundedness is a
                         running assumption of this type)

### Notes

This structure is optimized to evaluate the support function/vector with a large
sequence of directions that are close to each other. The strategy is to have an
index that can be used to warm-start the search for optimal values in the
support vector computation.

The default constructor assumes that the given list of edges is sorted.
It *does not perform* any sorting.
Use `addconstraint!` to iteratively add the edges in a sorted way.

- `HPolygonOpt(constraints::Vector{LinearConstraint{<:Real}}, [ind]::Int=1)`
  -- default constructor with optional index
"""
mutable struct HPolygonOpt{N<:Real} <: AbstractHPolygon{N}
    constraints::Vector{LinearConstraint{N}}
    ind::Int

    # default constructor that applies sorting of the given constraints
    function HPolygonOpt{N}(constraints::Vector{LinearConstraint{N}},
                            ind::Int=1;
                            sort_constraints::Bool=true,
                            check_boundedness::Bool=false,
                            prune::Bool=true) where {N<:Real}
        if sort_constraints
            sorted_constraints = Vector{LinearConstraint{N}}()
            sizehint!(sorted_constraints, length(constraints))
            for ci in constraints
                addconstraint!(sorted_constraints, ci; prune=prune)
            end
            P = new{N}(sorted_constraints, ind)
        else
            P = new{N}(constraints, ind)
        end
        @assert (!check_boundedness ||
                 isbounded(P, false)) "the polygon is not bounded"
        return P
    end
end

# convenience constructor without type parameter
HPolygonOpt(constraints::Vector{LinearConstraint{N}},
            ind::Int=1;
            sort_constraints::Bool=true,
            check_boundedness::Bool=false,
            prune::Bool=true) where {N<:Real} =
    HPolygonOpt{N}(constraints,
                   ind;
                   sort_constraints=sort_constraints,
                   check_boundedness=check_boundedness,
                   prune=prune)

# constructor with no constraints
HPolygonOpt{N}() where {N<:Real} = HPolygonOpt{N}(Vector{LinearConstraint{N}}())

# constructor with no constraints of type Float64
HPolygonOpt() = HPolygonOpt{Float64}()

# constructor from a simple H-representation
HPolygonOpt(A::AbstractMatrix{N},
            b::AbstractVector{N};
            sort_constraints::Bool=true,
            check_boundedness::Bool=false,
            prune::Bool=true) where {N<:Real} =
    HPolygonOpt(constraints_list(A, b); sort_constraints=sort_constraints,
                check_boundedness=check_boundedness, prune=prune)

# constructor from a simple H-representation with type parameter
HPolygonOpt{N}(A::AbstractMatrix{N},
               b::AbstractVector{N};
               sort_constraints::Bool=true,
               check_boundedness::Bool=false,
               prune::Bool=true) where {N<:Real} =
    HPolygonOpt(A, b; sort_constraints=sort_constraints,
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
