import Base.∈

export ZeroSet

"""
    ZeroSet <: LazySet

Type that represents the zero set, i.e., the set that only contains the origin.

### Fields

- `dim` -- the ambient dimension of this zero set
"""
struct ZeroSet <: LazySet
    dim::Int
end

"""
    dim(Z::ZeroSet)::Int

Return the ambient dimension of this zero set.

### Input

- `Z` -- a zero set, i.e., a set that only contains the origin

### Output

The ambient dimension of the zero set.
"""
function dim(Z::ZeroSet)::Int
    return Z.dim
end

"""
    σ(d, Z)

Return the support vector of a zero set.

### Input

- `Z` -- a zero set, i.e., a set that only contains the origin

### Output

The returned value is the origin since it is the only point that belongs to this
set.
"""
function σ(d::AbstractVector{N}, Z::ZeroSet)::Vector{N} where {N<:Real}
    return zeros(d)
end

"""
    ∈(x::AbstractVector, Z::ZeroSet)::Bool

Check whether a given point is contained in a zero set.

### Input

- `x` -- point/vector
- `Z` -- zero set

### Output

`true` iff ``x ∈ Z``.

### Examples

```jldoctest
julia> Z = ZeroSet(2);

julia> ∈([1.0, 0.0], Z)
false
julia> ∈([0.0, 0.0], Z)
true
```
"""
function ∈(x::AbstractVector{N}, Z::ZeroSet)::Bool where {N<:Real}
    @assert length(x) == dim(Z)

    zero_N = zero(N)
    return all(i -> x[i] == zero_N, eachindex(x))
end
