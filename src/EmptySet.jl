export EmptySet, ∅

"""
    EmptySet <: LazySet

Type that represents the empty set, i.e. the set with no elements.
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
"""
dim(S::EmptySet)::Int = -1

"""
    σ(d, ∅)

Return the support vector of an empty set.

### Input

- `∅` -- an empty set
"""
function σ(d::AbstractVector, S::EmptySet)
    error("the support vector of an empty set does not exist")
end
