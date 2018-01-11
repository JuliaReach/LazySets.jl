export Intersection, ∩

const EMPTY = 1::Int8
const NONEMPTY = 0::Int8
const UNKNOWN = -1::Int8

"""
    Intersection{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the intersection of two convex sets.

### Fields

- `X` -- convex set
- `Y` -- convex set
- `intersection` -- status if the intersection is empty or not (initially
                    unknown)
"""
struct Intersection{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2
    intersection::Int8

    # default constructor with dimension check
    Intersection{N, S1, S2}(X::S1, Y::S2) where
        {S1<:LazySet{N}, S2<:LazySet{N}} where {N<:Real} =
            dim(X) != dim(Y) ? throw(DimensionMismatch) : new(X, Y, UNKNOWN)
end
# type-less convenience constructor
Intersection(X::S1, Y::S2) where {S1<:LazySet{N}, S2<:LazySet{N}} where {N<:Real} =
    Intersection{N, S1, S2}(X, Y)

"""
    ∩

Alias for `Intersection`.
"""
function ∩(X, Y) = Intersection(X, Y)

"""
    Intersection(X, ∅)

Intersection of a set with the empty set from the right.

### Input

- `X` -- a convex set
- `∅` -- an empty set

### Output

An empty set, because the empty set is the absorbing element for the
intersection.
"""
Intersection(::LazySet, ∅::EmptySet) = ∅

"""
    Intersection(∅, X)

Intersection of a set with the empty set from the left.

### Input

- `∅` -- an empty set
- `X` -- a convex set

### Output

An empty set, because the empty set is the absorbing element for the
intersection.
"""
Intersection(∅::EmptySet, ::LazySet) = ∅

# special case: pure empty set intersection (we require the same numeric type)
(Intersection(∅::E, ::E)) where {E<:EmptySet} = ∅

"""
    dim(cap::Intersection)::Int

Return the dimension of an intersection of two convex sets.

### Input

- `cap` -- intersection of two convex sets

### Output

The ambient dimension of the intersection of two convex sets.
"""
function dim(cap::Intersection)::Int
    return dim(cap.X)
end

"""
    σ(d::AbstractVector{<:Real}, cap::Intersection)::AbstractVector{<:Real}

Return the support vector of an intersection of two convex sets in a given
direction.

### Input

- `d`   -- direction
- `cap` -- intersection of two convex sets

### Output

The support vector in the given direction.
If the direction has norm zero, the result TODO.
"""
function σ(d::AbstractVector{<:Real}, cap::Intersection)::AbstractVector{<:Real}
    if cap.intersection == EMPTY
        error("empty intersection, hence there is no support vector")
    end
    # TODO
    error("not implemented yet")
end

"""
    isempty(cap::Intersection)::Bool

Return if the intersection is empty or not.

### Input

- `cap` -- intersection of two convex sets

### Output

`true` iff the intersection is empty.
"""
function isempty(cap::Intersection)::Bool
    if cap.intersection != UNKNOWN
        return cap.intersection == EMPTY
    end
    # TODO
    error("not implemented yet")
end
