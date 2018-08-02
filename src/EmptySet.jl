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

"""
    norm(S::EmptySet, [p]::Real=Inf)

Return the norm of an empty set.
It is the norm of the enclosing ball (of the given ``p``-norm) of minimal volume
that is centered in the origin.

### Input

- `S` -- empty set
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function norm(S::EmptySet, p::Real=Inf)
    error("an empty set does not have a norm")
end

"""
    radius(S::EmptySet, [p]::Real=Inf)

Return the radius of an empty set.
It is the radius of the enclosing ball (of the given ``p``-norm) of minimal
volume with the same center.

### Input

- `S` -- empty set
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function radius(S::EmptySet, p::Real=Inf)
    error("an empty set does not have a radius")
end

"""
    diameter(S::EmptySet, [p]::Real=Inf)

Return the diameter of an empty set.
It is the maximum distance between any two elements of the set, or,
equivalently, the diameter of the enclosing ball (of the given ``p``-norm) of
minimal volume with the same center.

### Input

- `S` -- empty set
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function diameter(S::EmptySet, p::Real=Inf)
    error("an empty set does not have a diameter")
end
