import Base: isempty, ∈, ∩

export Intersection,
       IntersectionArray,
       array

"""
    Intersection{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the intersection of two convex sets.

### Fields

- `X` -- convex set
- `Y` -- convex set

### Examples

Create an expression, ``Z``, that lazily represents the intersection of two squares
``X`` and ``Y``:

```jldoctest lazy_intersection
julia> X, Y = BallInf([0,0.], 0.5), BallInf([1,0.], 0.65);

julia> Z = X ∩ Y;

julia> typeof(Z)
Intersection{Float64,BallInf{Float64},BallInf{Float64}}

julia> dim(Z)
2
```

We can check if the intersection is empty with `isempty`:

```jldoctest lazy_intersection
julia> isempty(Z)
false
```

Do not confuse `Intersection` with the concrete operation, that is computed with
the lowercase `intersection`:

```jldoctest lazy_intersection
julia> W = intersection(X, Y)
Hyperrectangle{Float64}([0.425, 0.0], [0.075, 0.5])
```
"""
struct Intersection{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    function Intersection{N, S1, S2}(X::S1, Y::S2) where
            {N<:Real, S1<:LazySet{N}, S2<:LazySet{N}}
        @assert dim(X) == dim(Y) "sets in an intersection must have the same " *
            "dimension"
        return new{N, S1, S2}(X, Y)
    end
end

# convenience constructor without type parameter
Intersection(X::S1, Y::S2) where {N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} =
    Intersection{N, S1, S2}(X, Y)

# EmptySet is the absorbing element for Intersection
@absorbing(Intersection, EmptySet)

"""
    ∩

Alias for `Intersection`.
"""
∩(X::LazySet, Y::LazySet) = Intersection(X, Y)


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
    σ(d::AbstractVector{N}, cap::Intersection{N}) where {N<:Real}

Return the support vector of an intersection of two convex sets in a given
direction.

### Input

- `d`   -- direction
- `cap` -- intersection of two convex sets

### Output

The support vector in the given direction.
"""
function σ(d::AbstractVector{N}, cap::Intersection{N}) where {N<:Real}
    # TODO document behavior if the direction has norm zero
    # TODO error message if the intersection is empty!
    # TODO implement
    error("not implemented yet")
end

"""
    ∈(x::AbstractVector{N}, cap::Intersection{N})::Bool where {N<:Real}

Check whether a given point is contained in an intersection of two convex sets.

### Input

- `x`   -- point/vector
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


# ================================
# intersection of an array of sets
# ================================

"""
    IntersectionArray{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the intersection of a finite number of convex sets.

### Fields

- `array` -- array of convex sets

### Notes

This type assumes that the dimensions of all elements match.

The `EmptySet` is the absorbing element for `IntersectionArray`.

Constructors:

- `IntersectionArray(array::Vector{<:LazySet})` -- default constructor

- `IntersectionArray([n]::Int=0, [N]::Type=Float64)`
  -- constructor for an empty sum with optional size hint and numeric type
"""
struct IntersectionArray{N<:Real, S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

@static if VERSION < v"0.7-"
    # convenience constructor without type parameter
    IntersectionArray(arr::Vector{S}) where {S<:LazySet{N}} where {N<:Real} =
        IntersectionArray{N, S}(arr)
end

# constructor for an empty sum with optional size hint and numeric type
function IntersectionArray(n::Int=0, N::Type=Float64)::IntersectionArray
    arr = Vector{LazySet{N}}()
    sizehint!(arr, n)
    return IntersectionArray(arr)
end

# EmptySet is the absorbing element for IntersectionArray
@absorbing(IntersectionArray, EmptySet)

# add functions connecting Intersection and IntersectionArray
@declare_array_version(Intersection, IntersectionArray)

"""
    array(ia::IntersectionArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}

Return the array of an intersection of a finite number of convex sets.

### Input

- `ia` -- intersection of a finite number of convex sets

### Output

The array of an intersection of a finite number of convex sets.
"""
function array(ia::IntersectionArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}
    return ia.array
end


# --- LazySet interface functions ---


"""
    dim(ia::IntersectionArray)::Int

Return the dimension of an intersection of a finite number of sets.

### Input

- `ia` -- intersection of a finite number of convex sets

### Output

The ambient dimension of the intersection of a finite number of sets.
"""
function dim(ia::IntersectionArray)::Int
    return length(ia.array) == 0 ? 0 : dim(ia.array[1])
end

"""
    σ(d::AbstractVector{<:Real}, ia::IntersectionArray)::Vector{<:Real}

Return the support vector of an intersection of a finite number of sets in a
given direction.

### Input

- `d`  -- direction
- `ia` -- intersection of a finite number of convex sets

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the individual sets.
"""
function σ(d::AbstractVector{<:Real}, ia::IntersectionArray)::Vector{<:Real}
    # TODO implement
    error("not implemented yet")
end

"""
    ∈(x::AbstractVector{N}, ia::IntersectionArray{N})::Bool where {N<:Real}

Check whether a given point is contained in an intersection of a finite number
of convex sets.

### Input

- `x`  -- point/vector
- `ia` -- intersection of a finite number of convex sets

### Output

`true` iff ``x ∈ ia``.
"""
function ∈(x::AbstractVector{N}, ia::IntersectionArray{N})::Bool where {N<:Real}
    for S in ia.array
        if x ∉ S
            return false
        end
    end
    return true
end
