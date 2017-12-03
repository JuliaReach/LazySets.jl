import Base.+

export MinkowskiSum,
       MinkowskiSumArray

"""
    MinkowskiSum{T1<:LazySet, T2<:LazySet} <: LazySet

Type that represents the Minkowski sum of two convex sets.

### Fields

- `X` -- a convex set
- `Y` -- a convex set
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
    +(X::LazySet, Y::LazySet)

Return the Minkowski sum of two convex sets.

### Input

- `X` -- convex set
- `Y` -- convex set

### Output

The Minkowski sum of the two convex sets.
"""
function +(X::LazySet, Y::LazySet)
    return MinkowskiSum(X, Y)
end

"""
    dim(ms::MinkowskiSum)

Return the dimension of a Minkowski sum.

### Input

- `ms` -- Minkowski sum

### Output

The ambient dimension of the Minkowski sum.
"""
function dim(ms::MinkowskiSum)
    return dim(ms.X)
end

"""
    σ(d::AbstractVector{<:Real}, ms::MinkowskiSum)::AbstractVector{<:Real}

Support vector of a Minkowski sum.

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

- `MinkowskiSumArray(n::Integer)` -- constructor for an empty sum with size hint
"""
struct MinkowskiSumArray{T<:LazySet} <: LazySet
    sfarray::Vector{T}
end
# constructor for an empty sum
MinkowskiSumArray() = MinkowskiSumArray{LazySet}(Vector{LazySet}(0))

function MinkowskiSumArray(n::Integer)::MinkowskiSumArray
    arr = Vector{LazySet}(0)
    sizehint!(arr, n)
    return MinkowskiSumArray(arr)
end

"""
    +(msa::MinkowskiSumArray, S::LazySet)::MinkowskiSumArray

Add a convex set to a Minkowski sum of a finite number of convex sets.

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

Add a convex set to a Minkowski sum of a finite number of convex sets.

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
    +(msa, vs)

Returns the original array because addition with a void set is a no-op.

### Input

- `msa` -- Minkowski sum array
- `vs` -- void set
"""
function +(msa::MinkowskiSumArray, ::VoidSet)::MinkowskiSumArray
    return msa
end

"""
    dim(msa::MinkowskiSumArray)

Return the dimension of a Minkowski sum of a finite number of sets.

### Input

- `msa` -- Minkowski sum array

### Output

The ambient dimension of the Minkowski sum of a finite number of sets.
"""
function dim(msa::MinkowskiSumArray)
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
