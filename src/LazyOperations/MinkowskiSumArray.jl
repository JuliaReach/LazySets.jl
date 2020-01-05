export MinkowskiSumArray,
       array

# ==================================
# Minkowski sum of an array of sets
# ==================================

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

isoperationtype(::Type{<:MinkowskiSumArray}) = true
isconvextype(::Type{MinkowskiSumArray{N, S}}) where {N, S} = isconvextype(S)

# constructor for an empty sum with optional size hint and numeric type
function MinkowskiSumArray(n::Int=0, N::Type=Float64)
   arr = Vector{LazySet{N}}()
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
   array(msa::MinkowskiSumArray{N, S}) where {N<:Real, S<:LazySet{N}}

Return the array of a Minkowski sum of a finite number of convex sets.

### Input

- `msa` -- Minkowski sum array

### Output

The array of a Minkowski sum of a finite number of convex sets.
"""
function array(msa::MinkowskiSumArray{N, S}) where {N<:Real, S<:LazySet{N}}
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
   σ(d::AbstractVector{N}, msa::MinkowskiSumArray{N}) where {N<:Real}

Return the support vector of a Minkowski sum of a finite number of sets in a
given direction.

### Input

- `d`   -- direction
- `msa` -- Minkowski sum array

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the summand sets.
"""
function σ(d::AbstractVector{N}, msa::MinkowskiSumArray{N}) where {N<:Real}
   return σ_helper(d, msa.array)
end

"""
   ρ(d::AbstractVector{N}, msa::MinkowskiSumArray{N}) where {N<:Real}

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
function ρ(d::AbstractVector{N}, msa::MinkowskiSumArray{N}) where {N<:Real}
   return sum([ρ(d, Xi) for Xi in msa.array])
end

"""
	isbounded(msa::MinkowskiSumArray)

Determine whether a Minkowski sum of a finite number of convex sets is bounded.

### Input

- `msa` -- Minkowski sum of a finite number of convex sets

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(msa::MinkowskiSumArray)
   return all(x -> isbounded(x), msa.array)
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
   return any(X -> isempty(X), array(msa))
end
