export MinkowskiSumArray,
       array

"""
   MinkowskiSumArray{N, S<:LazySet{N}} <: LazySet{N}

Type that represents the Minkowski sum of a finite number of sets.

### Fields

- `array` -- array of sets

### Notes

This type assumes that the dimensions of all elements match.

The `ZeroSet` is the neutral element and the `EmptySet` is the absorbing element
for `MinkowskiSumArray`.

The Minkowski sum preserves convexity: if the set arguments are convex, then
their Minkowski sum is convex as well.
"""
struct MinkowskiSumArray{N,S<:LazySet{N}} <: LazySet{N}
    array::Vector{S}
end

"""
    +(X::LazySet, Xs::LazySet...)
    +(Xs::Vector{<:LazySet})

Alias for the n-ary Minkowski sum.
"""
+(X::LazySet, Xs::LazySet...) = MinkowskiSumArray(vcat(X, Xs...))
+(X::LazySet) = X
+(Xs::Vector{<:LazySet}) = MinkowskiSumArray(Xs)

"""
    ⊕(X::LazySet, Xs::LazySet...)
    ⊕(Xs::Vector{<:LazySet})

Alias for the n-ary Minkowski sum.

### Notes

The function symbol can be typed via `\\oplus[TAB]`.
"""
⊕(X::LazySet, Xs::LazySet...) = +(X, Xs...)
⊕(Xs::Vector{<:LazySet}) = +(Xs)

isoperationtype(::Type{<:MinkowskiSumArray}) = true
isconvextype(::Type{MinkowskiSumArray{N,S}}) where {N,S} = isconvextype(S)

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
   array(msa::MinkowskiSumArray)

Return the array of a Minkowski sum of a finite number of sets.

### Input

- `msa` -- Minkowski sum of a finite number of sets

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

- `msa` -- Minkowski sum of a finite number of sets

### Output

The ambient dimension of the Minkowski sum of a finite number of sets, or `0` if
there is no set in the array.
"""
function dim(msa::MinkowskiSumArray)
    return length(msa.array) == 0 ? 0 : dim(msa.array[1])
end

"""
   σ(d::AbstractVector, msa::MinkowskiSumArray)

Return a support vector of a Minkowski sum of a finite number of sets in a given
direction.

### Input

- `d`   -- direction
- `msa` -- Minkowski sum of a finite number of sets

### Output

A support vector in the given direction.
If the direction has norm zero, the result depends on the summand sets.
"""
function σ(d::AbstractVector, msa::MinkowskiSumArray)
    return _σ_msum_array(d, msa.array)
end

@inline function _σ_msum_array(d::AbstractVector{N},
                               array::AbstractVector{<:LazySet}) where {N}
    return sum(σ(d, Xi) for Xi in array)
end

"""
   ρ(d::AbstractVector, msa::MinkowskiSumArray)

Evaluate the support function of a Minkowski sum of a finite number of sets in a
given direction.

### Input

- `d`   -- direction
- `msa` -- Minkowski sum of a finite number of sets

### Output

The evaluation of the support function in the given direction.

### Algorithm

The support function of the Minkowski sum of multiple sets evaluations to the
sum of the support-function evaluations of each set.
"""
function ρ(d::AbstractVector, msa::MinkowskiSumArray)
    return sum(ρ(d, Xi) for Xi in msa.array)
end

"""
    isbounded(msa::MinkowskiSumArray)

Check whether a Minkowski sum of a finite number of sets is bounded.

### Input

- `msa` -- Minkowski sum of a finite number of sets

### Output

`true` iff all wrapped sets are bounded.
"""
function isbounded(msa::MinkowskiSumArray)
    return all(isbounded, msa.array)
end

function isboundedtype(::Type{<:MinkowskiSumArray{N,S}}) where {N,S}
    return isboundedtype(S)
end

"""
   isempty(msa::MinkowskiSumArray)

Check whether a Minkowski sum of a finite number of sets is empty.

### Input

- `msa` -- Minkowski sum of a finite number of sets

### Output

`true` iff any of the wrapped sets is empty.
"""
function isempty(msa::MinkowskiSumArray)
    return any(isempty, array(msa))
end

"""
    center(msa::MinkowskiSumArray)

Return the center of a Minkowski sum of a finite number of centrally-symmetric
sets.

### Input

- `msa` -- Minkowski sum of a finite number of centrally-symmetric sets

### Output

The center of the set.
"""
function center(msa::MinkowskiSumArray)
    return sum(center(X) for X in msa)
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
