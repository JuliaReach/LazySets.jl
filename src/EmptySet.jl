import Base.∈

export EmptySet, ∅

"""
    EmptySet <: LazySet

Type that represents the empty set, i.e., the set with no elements.
"""
struct EmptySet <: LazySet end

"""
    ∅

An `EmptySet` instance.
"""
const ∅ = EmptySet()

"""
    dim(S::EmptySet)

Return the dimension of the empty set, which is -1 by convention.

### Input

- `S` -- an empty set

### Output

`-1` by convention.
"""
function dim(S::EmptySet)::Int
    return -1
end

"""
    σ(d, ∅)

Return the support vector of an empty set.

### Input

- `∅` -- an empty set

### Output

An error.
"""
function σ(d::AbstractVector, S::EmptySet)
    error("the support vector of an empty set does not exist")
end

"""
    ∈(x::AbstractVector, ∅::EmptySet)::Bool

Check whether a given point is contained in an empty set.

### Input

- `x` -- point/vector
- `∅` -- empty set

### Output

The output is always `false`.

### Examples

```jldoctest
julia> ∈([1.0, 0.0], ∅)
false
```
"""
function ∈(x::AbstractVector{<:Real}, Z::EmptySet)::Bool
    return false
end
