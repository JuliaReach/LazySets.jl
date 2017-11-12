export EmptySet

import Base:+,*

"""
    EmptySet <: LazySet

Type that represents the empty set, i.e. the set with no elements.
"""
struct EmptySet <: LazySet end

"""
    dim(∅::EmptySet)

Return the dimension of the empty set, which is -1 by convention.

### Input

- `∅` -- an empty set
"""
function dim(∅::EmptySet)
    return -1
end

"""
    σ(d, ∅)

Return the support vector of an empty set.

### Input

- `∅` -- an empty set
"""
function σ(d::AbstractVector{Float64}, ∅::EmptySet)::Vector{Float64}
    error("The support vector of an empty set does not exist")
end
