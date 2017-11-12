export ZeroSet

import Base:+,*

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
    dim(Z::ZeroSet)

Return the ambient dimension of this zero set.

### Input

- `Z` -- a zero set, i.e. a set which only contains the origin
"""
function dim(Z::ZeroSet)
    return Z.dim
end

"""
    σ(d, Z)

Return the support vector of a zero set.

### Input

- `Z` -- a zero set, i.e. a set which only contains the origin

### Output

The returned value is the origin since it's the only point that belongs to this
set.
"""
function σ(d::AbstractVector{Float64}, Z::ZeroSet)::Vector{Float64}
    return zeros(d)
end
