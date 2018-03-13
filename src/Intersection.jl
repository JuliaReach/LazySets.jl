import Base: isempty, ∈, ∩

export Intersection

"""
    Intersection{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the intersection of two convex sets.

### Fields

- `X` -- convex set
- `Y` -- convex set
"""
struct Intersection{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    Intersection{N, S1, S2}(X::S1, Y::S2) where
        {S1<:LazySet{N}, S2<:LazySet{N}} where {N<:Real} =
            dim(X) != dim(Y) ? throw(DimensionMismatch) : new(X, Y)
end
# type-less convenience constructor
Intersection(X::S1, Y::S2) where {S1<:LazySet{N}, S2<:LazySet{N}} where {N<:Real} =
    Intersection{N, S1, S2}(X, Y)

# EmptySet is the absorbing element for Intersection
@absorbing(Intersection, EmptySet)

"""
    ∩

Alias for `Intersection`.
"""
∩(X, Y) = Intersection(X, Y)


# --- LazySet interface functions ---


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
    σ(d::V, cap::Intersection)::V where {N<:Real, V<:AbstractVector{N}}

Return the support vector of an intersection of two convex sets in a given
direction.

### Input

- `d`   -- direction
- `cap` -- intersection of two convex sets

### Output

The support vector in the given direction.
"""
function σ(d::V, cap::Intersection)::V where {N<:Real, V<:AbstractVector{N}}
    # TODO document behavior if the direction has norm zero
    # TODO error message if the intersection is empty!
    # TODO implement
    error("not implemented yet")
end

"""
    ∈(x::AbstractVector{N}, cap::Intersection{N})::Bool where {N<:Real}

Check whether a given point is contained in a intersection of two convex sets.

### Input

- `x` -- point/vector
- `cap` -- intersection of two convex sets

### Output

`true` iff ``x ∈ cap``.
"""
function ∈(x::AbstractVector{N}, cap::Intersection{N})::Bool where {N<:Real}
    return (x ∈ cap.X) && (x ∈ cap.Y)
end


# --- Intersection functions ---


"""
    isempty(cap::Intersection)::Bool

Return if the intersection is empty or not.

### Input

- `cap` -- intersection of two convex sets

### Output

`true` iff the intersection is empty.
"""
function isempty(cap::Intersection)::Bool
    return is_intersection_empty(cap.X, cap.Y)
end
