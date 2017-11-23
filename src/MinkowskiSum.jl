import Base.+

export MinkowskiSum, MinkowskiSumArray

"""
    MinkowskiSum <: LazySet

Type that represents the Minkowski sum of two convex sets.

### Fields

- `X` -- a convex set
- `Y` -- a convex set
"""
struct MinkowskiSum{T1<:LazySet,T2<:LazySet} <: LazySet
    X::T1
    Y::T2

    # default constructor with dimension match check
    MinkowskiSum{T1,T2}(X::T1, Y::T2) where {T1<:LazySet,T2<:LazySet} =
        dim(X) != dim(Y) ? throw(DimensionMismatch) : new(X, Y)
end
# type-less convenience constructor
MinkowskiSum(X::T1, Y::T2) where {T1<:LazySet,T2<:LazySet} = MinkowskiSum{T1,T2}(X, Y)

"""
    X + Y

Convenience constructor for Minkowski sum.

### Input

- `X` -- a convex set
- `Y` -- another convex set

### Output

The symbolic Minkowski sum of ``X`` and ``Y``.
"""
function +(X::LazySet, Y::LazySet)
    return MinkowskiSum(X, Y)
end

"""
    dim(ms)

Ambient dimension of a Minkowski sum.

### Input

- `ms` -- Minkowski sum
"""
function dim(ms::MinkowskiSum)::Int64
    return dim(ms.X)
end

"""
    σ(d, ms)

Support vector of a Minkowski sum.

### Input

- `d`  -- vector
- `ms` -- Minkowski sum
"""
function σ(d::AbstractVector{<:Real}, ms::MinkowskiSum)::Vector{<:Real}
    return σ(d, ms.X) + σ(d, ms.Y)
end


# =================================
# Minkowski sum of an array of sets
# =================================
"""
    MinkowskiSumArray <: LazySet

Type that represents the Minkowski sum of a finite number of sets.

### Fields

- `sfarray` -- array of sets

### Notes

This type is optimized to be used on the left-hand side of additions only.
"""
struct MinkowskiSumArray{T<:LazySet} <: LazySet
    sfarray::Vector{T}

    MinkowskiSumArray{T}(sfarray::Vector{T}) where {T<:LazySet} = new(sfarray)
end
MinkowskiSumArray() = MinkowskiSumArray{LazySet}(Vector{LazySet}(0))
MinkowskiSumArray(sfarray::Vector{T}) where {T<:LazySet} = MinkowskiSumArray{T}(sfarray)

function MinkowskiSumArray(n::Int64)::MinkowskiSumArray
    arr = Vector{LazySet}(0)
    sizehint!(arr, n)
    return MinkowskiSumArray(arr)
end

"""
    +(msa, sf)

Adds the support function to the array.

### Input

- `msa` -- Minkowski sum array
- `sf` -- general support function

### Notes

This function is overridden for more specific types of `sf`.
"""
function +(msa::MinkowskiSumArray, sf::LazySet)::MinkowskiSumArray
    push!(msa.sfarray, sf)
    return msa
end

"""
    +(msa1, msa2)

Appends the elements of the second array to the first array.

### Input

- `msa1` -- first Minkowski sum array
- `msa2` -- second Minkowski sum array
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

"""
    dim(ms::MinkowskiSumArray)

Ambient dimension of the Minkowski sum of a finite number of sets.

### Input

- `ms` -- Minkowski sum array

### Notes

We do not double-check that the dimensions always match.
"""
function dim(ms::MinkowskiSumArray)::Int64
    return length(ms.sfarray) == 0 ? 0 : dim(ms.sfarray[1])
end

"""
    σ(d, ms)

Support vector of the Minkowski sum of a finite number of sets.

### Input

- `d` -- direction

- `ms` -- Minkowski sum array
"""
function σ(d::AbstractVector{<:Real}, ms::MinkowskiSumArray)::Vector{<:Real}
    svec = zeros(eltype(d), length(d))
    for sj in ms.sfarray
        svec += σ(d, sj)
    end
    return svec
end
