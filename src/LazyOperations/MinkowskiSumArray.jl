export MinkowskiSumArray,
       array

# ==================================
# Minkowski sum of an array of sets
# ==================================

"""
   MinkowskiSumArray{N, S<:ConvexSet{N}} <: ConvexSet{N}

Type that represents the Minkowski sum of a finite number of sets.

### Fields

- `array` -- array of sets

### Notes

This type assumes that the dimensions of all elements match.

The `ZeroSet` is the neutral element and the `EmptySet` is the absorbing element
for `MinkowskiSumArray`.

The Minkowski sum preserves convexity: if the set arguments are convex, then
their Minkowski sum is convex as well.

Constructors:

- `MinkowskiSumArray(array::Vector{<:ConvexSet})` -- default constructor

- `MinkowskiSumArray([n]::Int=0, [N]::Type=Float64)`
 -- constructor for an empty sum with optional size hint and numeric type
"""
struct MinkowskiSumArray{N, S<:ConvexSet{N}} <: ConvexSet{N}
   array::Vector{S}
end

isoperationtype(::Type{<:MinkowskiSumArray}) = true
isconvextype(::Type{MinkowskiSumArray{N, S}}) where {N, S} = isconvextype(S)

# constructor for an empty sum with optional size hint and numeric type
function MinkowskiSumArray(n::Int=0, N::Type=Float64)
   arr = Vector{ConvexSet{N}}()
   sizehint!(arr, n)
   return MinkowskiSumArray(arr)
end

# ZeroSet is the neutral element for MinkowskiSumArray
@neutral(MinkowskiSumArray, ZeroSet)

# EmptySet and Universe are the absorbing elements for MinkowskiSumArray
@absorbing(MinkowskiSumArray, EmptySet)
# @absorbing(MinkowskiSumArray, Universe)  # TODO problematic

# add functions connecting MinkowskiSum and MinkowskiSumArray
@declare_array_version(MinkowskiSum, MinkowskiSumArray)

"""
   array(msa::MinkowskiSumArray)

Return the array of a Minkowski sum of a finite number of sets.

### Input

- `msa` -- Minkowski sum array

### Output

The array of a Minkowski sum of a finite number of sets.
"""
function array(msa::MinkowskiSumArray)
   return msa.array
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
   return length(msa.array) == 0 ? 0 : dim(msa.array[1])
end

"""
   σ(d::AbstractVector, msa::MinkowskiSumArray)

Return the support vector of a Minkowski sum of a finite number of sets in a
given direction.

### Input

- `d`   -- direction
- `msa` -- Minkowski sum array

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the summand sets.
"""
function σ(d::AbstractVector, msa::MinkowskiSumArray)
   return _σ_msum_array(d, msa.array)
end

@inline function _σ_msum_array(d::AbstractVector{N},
                               array::AbstractVector{<:ConvexSet}) where {N}
    return sum(σ(d, Xi) for Xi in array)
end

"""
   ρ(d::AbstractVector, msa::MinkowskiSumArray)

Return the support function of a Minkowski sum array of a finite number of sets
in a given direction.

### Input

- `d`   -- direction
- `msa` -- Minkowski sum array

### Output

The support function in the given direction.

### Algorithm

The support function of the Minkowski sum of sets is the sum of the support
functions of each set.
"""
function ρ(d::AbstractVector, msa::MinkowskiSumArray)
   return sum(ρ(d, Xi) for Xi in msa.array)
end

"""
    isbounded(msa::MinkowskiSumArray)

Determine whether a Minkowski sum of a finite number of sets is bounded.

### Input

- `msa` -- Minkowski sum of a finite number of sets

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(msa::MinkowskiSumArray)
   return all(isbounded, msa.array)
end

function isboundedtype(::Type{<:MinkowskiSumArray{N, S}}) where {N, S}
    return isboundedtype(S)
end

"""
   isempty(msa::MinkowskiSumArray)

Return if a Minkowski sum array is empty or not.

### Input

- `msa` -- Minkowski sum array

### Output

`true` iff any of the wrapped sets are empty.
"""
function isempty(msa::MinkowskiSumArray)
   return any(isempty, array(msa))
end

"""
    center(msa::MinkowskiSumArray)

Return the center of a Minkowski sum array of centrally-symmetric sets.

### Input

- `msa` -- Minkowski sum array of centrally-symmetric sets

### Output

The center of the Minkowski sum array.
"""
function center(msa::MinkowskiSumArray)
    sum(center(X) for X in array(msa))
end

function concretize(msa::MinkowskiSumArray)
    a = array(msa)
    @assert !isempty(a) "an empty Minkowski sum is not allowed"
    X = msa
    @inbounds for (i, Y) in enumerate(a)
        if i == 1
            X = concretize(Y)
        else
            X = minkowski_sum(X, concretize(Y))
        end
    end
    return X
end
