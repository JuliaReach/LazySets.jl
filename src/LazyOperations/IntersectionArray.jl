export IntersectionArray,
       array

"""
    IntersectionArray{N, S<:LazySet{N}} <: LazySet{N}

Type that represents the intersection of a finite number of sets.

### Fields

- `array` -- array of sets

### Notes

This type assumes that the dimensions of all elements match.

The `EmptySet` is the absorbing element for `IntersectionArray`.

The intersection preserves convexity: if the set arguments are convex, then
their intersection is convex as well.

The convenience alias `∩` can be typed by `\\cap<tab>`.
"""
struct IntersectionArray{N,S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

"""
    Intersection!(X, Y)

Convenience function to compute the lazy intersection and modify
`IntersectionArray`s in-place.
"""
function Intersection! end

∩(X::LazySet, Xs::LazySet...) = IntersectionArray(vcat(X, Xs...))
∩(X::LazySet) = X
∩(Xs::Vector{<:LazySet}) = IntersectionArray(Xs)

isoperationtype(::Type{<:IntersectionArray}) = true
isconvextype(::Type{IntersectionArray{N,S}}) where {N,S} = isconvextype(S)

# constructor for an empty sum with optional size hint and numeric type
function IntersectionArray(n::Int=0, N::Type=Float64)
    arr = Vector{LazySet{N}}()
    sizehint!(arr, n)
    return IntersectionArray(arr)
end

# Universe is the neutral element for IntersectionArray
@neutral(IntersectionArray, Universe)

# EmptySet is the absorbing element for IntersectionArray
@absorbing(IntersectionArray, EmptySet)

# add functions connecting Intersection and IntersectionArray
@declare_array_version(Intersection, IntersectionArray)

"""
    array(ia::IntersectionArray)

Return the array of an intersection of a finite number of sets.

### Input

- `ia` -- intersection of a finite number of sets

### Output

The array of an intersection of a finite number of sets.
"""
function array(ia::IntersectionArray)
    return ia.array
end

"""
    dim(ia::IntersectionArray)

Return the dimension of an intersection of a finite number of sets.

### Input

- `ia` -- intersection of a finite number of sets

### Output

The ambient dimension of the intersection of a finite number of sets, or `0` if
there is no set in the array.
"""
function dim(ia::IntersectionArray)
    return length(ia.array) == 0 ? 0 : dim(ia.array[1])
end

"""
    σ(d::AbstractVector, ia::IntersectionArray)

Return a support vector of an intersection of a finite number of sets in a given
direction.

### Input

- `d`  -- direction
- `ia` -- intersection of a finite number of sets

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the individual sets.

### Algorithm

This implementation computes the concrete intersection, which can be expensive.
"""
@validate function σ(d::AbstractVector, ia::IntersectionArray)
    X = concretize(ia)
    return σ(d, X)
end

"""
    isbounded(ia::IntersectionArray)

Check whether an intersection of a finite number of sets is bounded.

### Input

- `ia` -- intersection of a finite number of sets

### Output

`true` iff the intersection is bounded.

### Algorithm

We first check if any of the wrapped sets is bounded.
Otherwise we check boundedness via
[`LazySets._isbounded_unit_dimensions`](@ref).
"""
function isbounded(ia::IntersectionArray)
    if any(isbounded, ia.array)
        return true
    end
    return _isbounded_unit_dimensions(ia)
end

function isboundedtype(::Type{<:IntersectionArray{N,S}}) where {N,S}
    return isboundedtype(S)
end

function ispolyhedraltype(::Type{<:IntersectionArray{N,S}}) where {N,S}
    return ispolyhedraltype(S)
end

function ispolyhedral(ia::IntersectionArray)
    ispolyhedraltype(typeof(ia)) && return true

    return all(ispolyhedral, array(ia))
end

"""
    ∈(x::AbstractVector, ia::IntersectionArray)

Check whether a given point is contained in an intersection of a finite number
of sets.

### Input

- `x`  -- point/vector
- `ia` -- intersection of a finite number of sets

### Output

`true` iff ``x ∈ ia``.

### Algorithm

A point ``x`` is in the intersection iff it is in each set.
"""
@validate function ∈(x::AbstractVector, ia::IntersectionArray)
    for S in ia.array
        if x ∉ S
            return false
        end
    end
    return true
end

"""
    constraints_list(ia::IntersectionArray)

Return the list of constraints of an intersection of a finite number of
(polyhedral) sets.

### Input

- `ia` -- intersection of a finite number of (polyhedral) sets

### Output

The list of constraints of the intersection.

### Notes

We assume that the underlying sets are polyhedral, i.e., offer a method
`constraints_list`.

### Algorithm

We create the polyhedron from the `constraints_list`s of the sets and remove
redundant constraints.
"""
@validate function constraints_list(ia::IntersectionArray)
    N = eltype(ia)
    constraints = Vector{HalfSpace{N,Vector{N}}}() # TODO: use vector type of ia
    for X in ia
        clist_X = _normal_Vector(X)
        if clist_X isa HalfSpace
            push!(constraints, clist_X)
        else
            append!(constraints, clist_X)
        end
    end
    remove_redundant_constraints!(constraints)
    return constraints
end

function concretize(ia::IntersectionArray)
    return _concretize_lazy_array(ia)
end

@validate function translate(ia::IntersectionArray, x::AbstractVector)
    return IntersectionArray([translate(X, x) for X in array(ia)])
end
