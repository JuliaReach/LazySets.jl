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

    MinkowskiSum{T1,T2}(X::T1, Y::T2) where {T1<:LazySet,T2<:LazySet} = dim(X) != dim(Y) ? throw(DimensionMismatch) : new(X, Y)
end
MinkowskiSum(X::T1, Y::T2) where {T1<:LazySet,T2<:LazySet} = MinkowskiSum{T1,T2}(X, Y)

function +(X::LazySet, Y::LazySet)
    return MinkowskiSum(X, Y)
end

# ambient dimension of the Minkowski sum of two sets
function dim(ms::MinkowskiSum)::Int64
    return dim(ms.X)
end

# support vector of the Minkowski sum of two sets
function σ(d::AbstractVector{Float64}, ms::MinkowskiSum)::Vector{Float64}
    return σ(d, ms.X) + σ(d, ms.Y)
end


# =====================
# Minkowski sum of sets
# =====================
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
    σ(d::Vector{Float64}, ms::MinkowskiSumArray)

Support vector of the Minkowski sum of a finite number of sets.

### Input

- `d` -- direction

- `ms` -- Minkowski sum array
"""
function σ(d::AbstractVector{Float64}, ms::MinkowskiSumArray)::Vector{Float64}
    svec = zeros(length(d))
    for sj in ms.sfarray
        svec += σ(d, sj)
    end
    return svec
end
