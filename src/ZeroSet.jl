import Base.∈

export ZeroSet,
       linear_map

"""
    ZeroSet{N<:Real} <: AbstractSingleton{N}

Type that represents the zero set, i.e., the set that only contains the origin.

### Fields

- `dim` -- the ambient dimension of this zero set
"""
struct ZeroSet{N<:Real} <: AbstractSingleton{N}
    dim::Int
end
# type-less convenience constructor
ZeroSet(dim::Int) = ZeroSet{Float64}(dim)


# --- AbstractSingleton interface functions ---


"""
    element(S::ZeroSet{N})::Vector{N} where {N<:Real}

Return the element of a zero set.

### Input

- `S` -- zero set

### Output

The element of the zero set, i.e., a zero vector.
"""
function element(S::ZeroSet{N})::Vector{N} where {N<:Real}
    return zeros(N, S.dim)
end

"""
    element(S::ZeroSet{N}, ::Int)::N where {N<:Real}

Return the i-th entry of the element of a zero set.

### Input

- `S` -- zero set
- `i` -- dimension

### Output

The i-th entry of the element of the zero set, i.e., 0.
"""
function element(S::ZeroSet{N}, ::Int)::N where {N<:Real}
    return zero(N)
end


# --- LazySet interface functions ---


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
    σ(d::AbstractVector{N}, Z::ZeroSet)::Vector{N} where {N<:Real}

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
    ∈(x::AbstractVector{N}, Z::ZeroSet{N})::Bool where {N<:Real}

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
function ∈(x::AbstractVector{N}, Z::ZeroSet{N})::Bool where {N<:Real}
    @assert length(x) == dim(Z)

    zero_N = zero(N)
    return all(i -> x[i] == zero_N, eachindex(x))
end

"""
    linear_map(M::AbstractMatrix, Z::ZeroSet{N}) where {N<:Real}

Concrete linear map of a zero set.

### Input

- `M` -- matrix
- `Z` -- zero set

### Output

The zero set whose dimension matches the output dimension of the given matrix.
"""
function linear_map(M::AbstractMatrix, Z::ZeroSet{N}) where {N<:Real}
    m, n = size(M)
    return ZeroSet(m)
end
