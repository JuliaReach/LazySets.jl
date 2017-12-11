import Base.+

export MinkowskiSum, ⊕,
       MinkowskiSumArray

"""
    MinkowskiSum{T1<:LazySet, T2<:LazySet} <: LazySet

Type that represents the Minkowski sum of two convex sets.

### Fields

- `X` -- first convex set
- `Y` -- second convex set
"""
struct MinkowskiSum{T1<:LazySet, T2<:LazySet} <: LazySet
    X::T1
    Y::T2

    # default constructor with dimension match check
    MinkowskiSum{T1,T2}(X::T1, Y::T2) where {T1<:LazySet,T2<:LazySet} =
        dim(X) != dim(Y) ? throw(DimensionMismatch) : new(X, Y)
end
# type-less convenience constructor
MinkowskiSum(X::T1, Y::T2) where {T1<:LazySet, T2<:LazySet} =
    MinkowskiSum{T1,T2}(X, Y)

"""
    X + Y

Convenience constructor for Minkowski sum.

### Input

- `X` -- a convex set
- `Y` -- another convex set

### Output

The symbolic Minkowski sum of ``X`` and ``Y``.
"""
+(X::LazySet, Y::LazySet) = MinkowskiSum(X, Y)

⊕(X::LazySet, Y::LazySet) = +(X, Y)

"""
    X + ∅

Right Minkowski sum of a set by an empty set.

### Input

- `X` -- a convex set
- `∅` -- an empty set

### Output

An empty set, because the empty set is the absorbing element for the
Minkowski sum.
"""
+(X::LazySet, ∅::EmptySet) = ∅

"""
    ∅ + X

Left Minkowski sum of a set by an empty set.

### Input

- `∅` -- an empty set
- `X` -- a convex set

### Output

An empty set, because the empty set is the absorbing element for the
Minkowski sum.
"""
+(∅::EmptySet, X::LazySet) = ∅

+(X::LazySet, ::ZeroSet) = X

+(::ZeroSet, X::LazySet) = X

"""
    dim(ms::MinkowskiSum)::Int

Return the dimension of a Minkowski sum.

### Input

- `ms` -- Minkowski sum

### Output

The ambient dimension of the Minkowski sum.
"""
dim(ms::MinkowskiSum)::Int = dim(ms.X)

"""
    σ(d::AbstractVector{<:Real}, ms::MinkowskiSum)::AbstractVector{<:Real}

Return the support vector of a Minkowski sum.

### Input

- `d`  -- direction
- `ms` -- Minkowski sum

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the summand sets.
"""
function σ(d::AbstractVector{<:Real}, ms::MinkowskiSum)::AbstractVector{<:Real}
    return σ(d, ms.X) + σ(d, ms.Y)
end


# =================================
# Minkowski sum of an array of sets
# =================================
"""
    MinkowskiSumArray{T<:LazySet} <: LazySet

Type that represents the Minkowski sum of a finite number of convex sets.

### Fields

- `sfarray` -- array of convex sets

### Notes

This type assumes that the dimensions of all elements match.

- `MinkowskiSumArray(sfarray::Vector{<:LazySet})` -- default constructor

- `MinkowskiSumArray()` -- constructor for an empty sum

- `MinkowskiSumArray(n::Int)` -- constructor for an empty sum with size hint
"""
struct MinkowskiSumArray{T<:LazySet} <: LazySet
    sfarray::Vector{T}
end
# constructor for an empty sum
MinkowskiSumArray() = MinkowskiSumArray{LazySet}(Vector{LazySet}(0))
# constructor for an empty sum with size hint
function MinkowskiSumArray(n::Int)::MinkowskiSumArray
    arr = Vector{LazySet}(0)
    sizehint!(arr, n)
    return MinkowskiSumArray(arr)
end

"""
    +(msa::MinkowskiSumArray, S::LazySet)::MinkowskiSumArray

Add a convex set to a Minkowski sum of a finite number of convex sets from the
right.

### Input

- `msa` -- Minkowski sum array (is modified)
- `S`   -- convex set

### Output

The modified Minkowski sum of a finite number of convex sets.
"""
function +(msa::MinkowskiSumArray, S::LazySet)::MinkowskiSumArray
    push!(msa.sfarray, S)
    return msa
end

"""
    +(S::LazySet, msa::MinkowskiSumArray)::MinkowskiSumArray

Add a convex set to a Minkowski sum of a finite number of convex sets from the
left.

### Input

- `S`   -- convex set
- `msa` -- Minkowski sum array (is modified)

### Output

The modified Minkowski sum of a finite number of convex sets.
"""
function +(S::LazySet, msa::MinkowskiSumArray)::MinkowskiSumArray
    push!(msa.sfarray, S)
    return msa
end

"""
    +(msa1::MinkowskiSumArray, msa2::MinkowskiSumArray)::MinkowskiSumArray

Add the elements of a finite Minkowski sum of convex sets to another finite
Minkowski sum.

### Input

- `msa1` -- first Minkowski sum array (is modified)
- `msa2` -- second Minkowski sum array

### Output

The modified first Minkowski sum of a finite number of convex sets.
"""
function +(msa1::MinkowskiSumArray, msa2::MinkowskiSumArray)::MinkowskiSumArray
    append!(msa1.sfarray, msa2.sfarray)
    return msa1
end

"""
    +(msa, Z)

Returns the original array because addition with an empty set is a no-op.

### Input

- `msa` -- Minkowski sum array
- `Z`  -- a Zero set
"""
function +(msa::MinkowskiSumArray, Z::ZeroSet)::MinkowskiSumArray
    return msa
end

function +(Z::ZeroSet, msa::MinkowskiSumArray)::MinkowskiSumArray
    return msa
end

"""
    +(msa::MinkowskiSumArray, ∅::EmptySet)

Right Minkowski sum of a set by an empty set.

### Input

- `msa` -- Minkowski sum array
- `∅`   -- an empty set

### Output

An empty set, because the empty set is the absorbing element for the
Minkowski sum.
"""
+(msa::MinkowskiSumArray, ∅::EmptySet) = ∅

"""
    +(∅::EmptySet, msa::MinkowskiSumArray)

Left Minkowski sum of a set by an empty set.

### Input

- `∅` -- an empty set
- `msa` -- Minkowski sum array

### Output

An empty set, because the empty set is the absorbing element for the
Minkowski sum.
"""
+(∅::EmptySet, msa::MinkowskiSumArray) = ∅

"""
    dim(msa::MinkowskiSumArray)::Int

Return the dimension of a Minkowski sum of a finite number of sets.

### Input

- `msa` -- Minkowski sum array

### Output

The ambient dimension of the Minkowski sum of a finite number of sets.
"""
function dim(msa::MinkowskiSumArray)::Int
    return length(msa.sfarray) == 0 ? 0 : dim(msa.sfarray[1])
end

"""
    σ(d::AbstractVector{<:Real}, msa::MinkowskiSumArray)::Vector{<:Real}

Return the support vector of a Minkowski sum of a finite number of sets in a
given direction.

### Input

- `d`   -- direction

- `msa` -- Minkowski sum array

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the summand sets.
"""
function σ(d::AbstractVector{<:Real}, msa::MinkowskiSumArray)::Vector{<:Real}
    svec = zeros(eltype(d), length(d))
    for sj in msa.sfarray
        svec += σ(d, sj)
    end
    return svec
end
