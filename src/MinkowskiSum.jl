import Base.+

export MinkowskiSum, ⊕,
       MinkowskiSumArray,
       array

"""
    MinkowskiSum{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the Minkowski sum of two convex sets.

### Fields

- `X` -- first convex set
- `Y` -- second convex set

### Notes

The `ZeroSet` is the neutral element and the `EmptySet` is the absorbing element
for `MinkowskiSum`.
"""
struct MinkowskiSum{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension match check
    MinkowskiSum{N, S1, S2}(X::S1, Y::S2) where
        {S1<:LazySet{N}, S2<:LazySet{N}} where {N<:Real} =
            dim(X) != dim(Y) ? throw(DimensionMismatch) : new(X, Y)
end
# type-less convenience constructor
MinkowskiSum(X::S1, Y::S2) where {S1<:LazySet{N}, S2<:LazySet{N}} where {N<:Real} =
    MinkowskiSum{N, S1, S2}(X, Y)

# ZeroSet is the neutral element and EmptySet is the absorbing element for
# MinkowskiSum
@commutative_neutral_absorbing(MinkowskiSum, ZeroSet, EmptySet)

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

"""
    ⊕(X::LazySet, Y::LazySet)

Unicode alias constructor `\oplus` for the Minkowski sum operator `+(X, Y)`.
"""
⊕(X::LazySet, Y::LazySet) = +(X, Y)

"""
    dim(ms::MinkowskiSum)::Int

Return the dimension of a Minkowski sum.

### Input

- `ms` -- Minkowski sum

### Output

The ambient dimension of the Minkowski sum.
"""
function dim(ms::MinkowskiSum)::Int
    return dim(ms.X)
end

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
    MinkowskiSumArray{N<:Real, S<:LazySet{N}} <: LazySet{N}

Type that represents the Minkowski sum of a finite number of convex sets.

### Fields

- `array` -- array of convex sets

### Notes

This type assumes that the dimensions of all elements match.

The `ZeroSet` is the neutral element and the `EmptySet` is the absorbing element
for `MinkowskiSumArray`.

Constructors:

- `MinkowskiSumArray(array::Vector{<:LazySet})` -- default constructor

- `MinkowskiSumArray([n]::Int=0, [N]::Type=Float64)`
  -- constructor for an empty sum with optional size hint and numeric type
"""
struct MinkowskiSumArray{N<:Real, S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

# type-less convenience constructor
MinkowskiSumArray(arr::Vector{S}) where {S<:LazySet{N}} where {N<:Real} =
    MinkowskiSumArray{N, S}(arr)

# constructor for an empty sum with optional size hint and numeric type
function MinkowskiSumArray(n::Int=0, N::Type=Float64)::MinkowskiSumArray
    arr = Vector{LazySet{N}}(0)
    sizehint!(arr, n)
    return MinkowskiSumArray(arr)
end

"""
    MinkowskiSumArray(msa::MinkowskiSumArray, S::LazySet)::MinkowskiSumArray

Add a convex set to a Minkowski sum of a finite number of convex sets from the
right.

### Input

- `msa` -- Minkowski sum array (is modified)
- `S`   -- convex set

### Output

The modified Minkowski sum of a finite number of convex sets.
"""
function MinkowskiSumArray(msa::MinkowskiSumArray,
                           S::LazySet)::MinkowskiSumArray
    push!(msa.array, S)
    return msa
end

"""
    MinkowskiSumArray(S::LazySet, msa::MinkowskiSumArray)::MinkowskiSumArray

Add a convex set to a Minkowski sum of a finite number of convex sets from the
left.

### Input

- `S`   -- convex set
- `msa` -- Minkowski sum array (is modified)

### Output

The modified Minkowski sum of a finite number of convex sets.
"""
function MinkowskiSumArray(S::LazySet,
                           msa::MinkowskiSumArray)::MinkowskiSumArray
    return MinkowskiSumArray(msa, S)
end

"""
    MinkowskiSumArray(msa1::MinkowskiSumArray,
                      msa2::MinkowskiSumArray)::MinkowskiSumArray

Add the elements of a finite Minkowski sum of convex sets to another finite
Minkowski sum.

### Input

- `msa1` -- first Minkowski sum array (is modified)
- `msa2` -- second Minkowski sum array

### Output

The modified first Minkowski sum of a finite number of convex sets.
"""
function MinkowskiSumArray(msa1::MinkowskiSumArray,
                           msa2::MinkowskiSumArray)::MinkowskiSumArray
    append!(msa1.array, msa2.array)
    return msa1
end

# ZeroSet is the neutral element and EmptySet is the absorbing element for
# MinkowskiSumArray
@commutative_neutral_absorbing(MinkowskiSumArray, ZeroSet, EmptySet)

"""
    array(msa::MinkowskiSumArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}

Return the array of a Minkowski sum of a finite number of convex sets.

### Input

- `msa` -- Minkowski sum array

### Output

The array of a Minkowski sum of a finite number of convex sets.
"""
function array(msa::MinkowskiSumArray{N, S})::Vector{S} where {N<:Real, S<:LazySet{N}}
    return msa.array
end

"""
    dim(msa::MinkowskiSumArray)::Int

Return the dimension of a Minkowski sum of a finite number of sets.

### Input

- `msa` -- Minkowski sum array

### Output

The ambient dimension of the Minkowski sum of a finite number of sets.
"""
function dim(msa::MinkowskiSumArray)::Int
    return length(msa.array) == 0 ? 0 : dim(msa.array[1])
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
    for sj in msa.array
        svec += σ(d, sj)
    end
    return svec
end
