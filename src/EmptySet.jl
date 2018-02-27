import Base.∈

export EmptySet, ∅,
       an_element

"""
    EmptySet{N<:Real} <: LazySet{N}

Type that represents the empty set, i.e., the set with no elements.
"""
struct EmptySet{N<:Real} <: LazySet{N} end

# type-less convenience constructor
EmptySet() = EmptySet{Float64}()

"""
    ∅

An `EmptySet` instance of type `Float64`.
"""
const ∅ = EmptySet{Float64}()

"""
    dim(∅::EmptySet)

Return the dimension of the empty set, which is -1 by convention.

### Input

- `∅` -- an empty set

### Output

`-1` by convention.
"""
function dim(∅::EmptySet)::Int
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
function σ(d::AbstractVector, ∅::EmptySet)
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
function ∈(x::AbstractVector{<:Real}, ∅::EmptySet)::Bool
    return false
end

"""
    an_element(∅::EmptySet)

Return some element of an empty set.

### Input

- `∅` -- empty set

### Output

An error.
"""
function an_element(∅::EmptySet)
    error("an empty set does not have any element")
end
