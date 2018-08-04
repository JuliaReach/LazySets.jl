import Base.<=

export HPolygonOpt

"""
    HPolygonOpt{N<:Real} <: AbstractHPolygon{N}

Type that represents a convex polygon in constraint representation whose edges
are sorted in counter-clockwise fashion with respect to their normal directions.
This is a refined version of `HPolygon`.

### Fields

- `constraints` -- list of linear constraints
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

- `HPolygonOpt(constraints::Vector{LinearConstraint{<:Real}}, [ind]::Int)`
  -- default constructor with optional index
- `HPolygonOpt(S::LazySet)` -- constructor from another set
"""
mutable struct HPolygonOpt{N<:Real} <: AbstractHPolygon{N}
    constraints::Vector{LinearConstraint{N}}
    ind::Int

    # default constructor
    HPolygonOpt{N}(constraints::Vector{LinearConstraint{N}},
                   ind::Int=1) where {N<:Real} =
        new{N}(constraints, ind)
end
# type-less convenience constructor with optional index
HPolygonOpt(constraints::Vector{LinearConstraint{N}},
            ind::Int=1) where {N<:Real} =
    HPolygonOpt{N}(constraints, ind)

# constructor for an HPolygon with no constraints
HPolygonOpt{N}() where {N<:Real} = HPolygonOpt{N}(Vector{LinearConstraint{N}}(0), 1)

# constructor for an HPolygon with no constraints of type Float64
HPolygonOpt() = HPolygonOpt{Float64}()

# conversion constructor
HPolygonOpt(S::LazySet) = convert(HPolygonOpt, S)



# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, P::HPolygonOpt)::Vector{N} where {N<:Real}

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
function σ(d::AbstractVector{N}, P::HPolygonOpt)::Vector{N} where {N<:Real}
    n = length(P.constraints)
    @assert n > 0 "the polygon has no constraints"
    if (d <= P.constraints[P.ind].a)
        k = P.ind-1
        while (k >= 1 && d <= P.constraints[k].a)
            k -= 1
        end
        if (k == 0)
            P.ind = n
            return element(intersection(Line(P.constraints[n]),
                                        Line(P.constraints[1])))
        else
            P.ind = k
            return element(intersection(Line(P.constraints[k]),
                                        Line(P.constraints[k+1])))
        end
    else
        k = P.ind+1
        while (k <= n && P.constraints[k].a <= d)
            k += 1
        end
        if (k == n+1)
            P.ind = n
            return element(intersection(Line(P.constraints[n]),
                                        Line(P.constraints[1])))
        else
            P.ind = k-1
            return element(intersection(Line(P.constraints[k-1]),
                                        Line(P.constraints[k])))
        end
    end
end
