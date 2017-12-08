export ZeroSet

"""
    ZeroSet <: LazySet

Type that represents the zero set, i.e. the set which only contains the origin.

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

- `Z` -- a zero set, i.e. a set which only contains the origin
"""
dim(Z::ZeroSet)::Int = Z.dim

"""
    σ(d, Z)

Return the support vector of a zero set.

### Input

- `Z` -- a zero set, i.e. a set which only contains the origin

### Output

The returned value is the origin since it's the only point that belongs to this
set.
"""
function σ(d::AbstractVector{<:Real}, Z::ZeroSet)::Vector{<:Real}
    return zeros(d)
end
