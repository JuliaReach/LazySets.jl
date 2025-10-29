"""
    HPolygon{N, VN<:AbstractVector{N}} <: AbstractHPolygon{N}

Type that represents a convex polygon in constraint representation whose edges
are sorted in counter-clockwise fashion with respect to their normal directions.

### Fields

- `constraints` -- list of linear constraints, sorted by the normal direction in
                   counter-clockwise fashion

### Notes

Further constructor arguments:

- `sort_constraints`  -- (optional, default: `true`) flag for sorting the
                         constraints (being sorted is a running assumption of
                         this type)
- `check_boundedness` -- (optional, default: `false`) flag for checking if the
                         constraints make the polygon bounded; (boundedness is a
                         running assumption of this type)
- `prune`             -- (optional, default: `true`) flag for removing redundant
                         constraints

The option `sort_constraints` can be used to deactivate automatic sorting of
constraints in counter-clockwise fashion, which is an invariant of this type.
Alternatively, one can construct an `HPolygon` with empty constraints list,
which can then be filled iteratively using `addconstraint!`.

Similarly, the option `prune` can be used to deactivate automatic pruning of
redundant constraints.

Another type assumption is that the polygon is bounded.
The option `check_boundedness` can be used to assert this.
This option is deactivated by default because we explicitly want to allow the
iterative addition of the constraints, and hence one has to initially construct
an empty list of constraints (which represents an unbounded set).
The user has to make sure that the `HPolygon` is not used before the constraints
actually describe a bounded set.
"""
struct HPolygon{N,VN<:AbstractVector{N}} <: AbstractHPolygon{N}
    constraints::Vector{HalfSpace{N,VN}}

    # default constructor that applies sorting of the given constraints and
    # (checks for and) removes redundant constraints
    function HPolygon(constraints::Vector{HalfSpace{N,VN}};
                      sort_constraints::Bool=true,
                      check_boundedness::Bool=false,
                      prune::Bool=true) where {N,VN<:AbstractVector{N}}
        if sort_constraints
            sorted_constraints = Vector{HalfSpace{N,VN}}()
            sizehint!(sorted_constraints, length(constraints))
            for ci in constraints
                addconstraint!(sorted_constraints, ci; prune=prune)
            end
            P = new{N,VN}(sorted_constraints)
        else
            P = new{N,VN}(constraints)
        end
        @assert (!check_boundedness || isbounded(P, false)) "the polygon is not bounded"
        return P
    end
end

# constructor with no constraints
function HPolygon{N,VN}() where {N,VN<:AbstractVector{N}}
    return HPolygon(Vector{HalfSpace{N,VN}}())
end

# constructor with no constraints and given numeric type
function HPolygon{N}() where {N}
    return HPolygon(Vector{HalfSpace{N,Vector{N}}}())
end

# constructor without explicit numeric type, defaults to Float64
function HPolygon()
    return HPolygon{Float64}()
end

# constructor with constraints of mixed type
function HPolygon(constraints::Vector{<:HalfSpace})
    return HPolygon(_normal_Vector(constraints))
end

# constructor from a simple constraint representation
function HPolygon(A::AbstractMatrix, b::AbstractVector; sort_constraints::Bool=true,
                  check_boundedness::Bool=false, prune::Bool=true)
    return HPolygon(constraints_list(A, b); sort_constraints=sort_constraints,
                    check_boundedness=check_boundedness, prune=prune)
end
