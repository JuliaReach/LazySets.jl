export IntersectionArray,
       array

# ================================
# Intersection of an array of sets
# ================================

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

Constructors:

- `IntersectionArray(array::Vector{<:LazySet})` -- default constructor

- `IntersectionArray([n]::Int=0, [N]::Type=Float64)`
 -- constructor for an empty sum with optional size hint and numeric type
"""
struct IntersectionArray{N, S<:LazySet{N}} <: LazySet{N}
   array::Vector{S}
end

isoperationtype(::Type{<:IntersectionArray}) = true
isconvextype(::Type{IntersectionArray{N, S}}) where {N, S} = isconvextype(S)

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


# --- LazySet interface functions ---


"""
    dim(ia::IntersectionArray)

Return the dimension of an intersection of a finite number of sets.

### Input

- `ia` -- intersection of a finite number of sets

### Output

The ambient dimension of the intersection of a finite number of sets.
"""
function dim(ia::IntersectionArray)
   return length(ia.array) == 0 ? 0 : dim(ia.array[1])
end

"""
    σ(d::AbstractVector, ia::IntersectionArray)

Return the support vector of an intersection of a finite number of sets in a
given direction.

### Input

- `d`  -- direction
- `ia` -- intersection of a finite number of sets

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the individual sets.
"""
function σ(d::AbstractVector, ia::IntersectionArray)
    X = concretize(ia)
    return σ(d, X)
end

"""
    isbounded(ia::IntersectionArray)

Determine whether an intersection of a finite number of sets is bounded.

### Input

- `ia` -- intersection of a finite number of sets

### Output

`true` iff the intersection is bounded.

### Algorithm

We first check if any of the wrapped sets is bounded.
Otherwise, we check boundedness via [`LazySets._isbounded_unit_dimensions`](@ref).
"""
function isbounded(ia::IntersectionArray)
   if any(isbounded, ia.array)
       return true
   end
   return _isbounded_unit_dimensions(ia)
end

function isboundedtype(::Type{<:IntersectionArray{N, S}}) where {N, S}
    return isboundedtype(S)
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
"""
function ∈(x::AbstractVector, ia::IntersectionArray)
   for S in ia.array
       if x ∉ S
           return false
       end
   end
   return true
end

"""
    constraints_list(ia::IntersectionArray{N}) where {N}

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
function constraints_list(ia::IntersectionArray{N}) where {N}
   constraints = Vector{LinearConstraint{N, Vector{N}}}() # TODO: use vector type of ia
   for X in array(ia)
       clist_X = _normal_Vector(X)
       append!(constraints, clist_X)
   end
   remove_redundant_constraints!(constraints)
   return constraints
end

function concretize(ia::IntersectionArray)
    a = array(ia)
    @assert !isempty(a) "an empty intersection is not allowed"
    X = ia
    @inbounds for (i, Y) in enumerate(a)
        if i == 1
            X = concretize(Y)
        else
            X = intersection(X, concretize(Y))
        end
    end
    return X
end
