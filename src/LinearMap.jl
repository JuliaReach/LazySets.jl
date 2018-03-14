import Base: *, ∈

export LinearMap

"""
    LinearMap{NM, N} <: LazySet{N}

Type that represents a linear transformation ``M⋅S`` of a convex set ``S``.

### Fields

- `M` -- matrix/linear map
- `X` -- convex set

### Notes

This type is parametric in the elements of the linear map, `NM`, and independently
on the type of elements of the target set (`N`). Typically `NM = N` but it is
not necessarily always the case e.g. if `NM` is an interval that holds numbers of
type `N`, where `N` is a floating point number type such as `Float64`.
"""
struct LinearMap{NM, N} <: LazySet{N}
    M::AbstractMatrix{NM}
    X::LazySet{N}

    # default constructor with dimension match check
    function LinearMap{NM, N}(M::AbstractMatrix{NM}, X::LazySet{N}) where {NM, N}
        @assert dim(X) == size(M, 2) "a linear map of size $(size(M)) cannot be applied to a set of dimension $(dim(X))"
        return new(M, X)
    end
end

# type-less convenience constructor
LinearMap(M::AbstractMatrix{NM}, X::LazySet{N}) where {NM, N} = LinearMap{NM, N}(M, X)

# constructor from a linear map: perform the matrix multiplication immediately
LinearMap(M::AbstractMatrix{NM},
          map::LinearMap{NM, N}) where {NM, N} = LinearMap{NM, N}(M * map.M, map.X)

"""
```
    *(M::AbstractMatrix, X::LazySet)
```

Return the linear map of a convex set.

### Input

- `M` -- matrix/linear map
- `X` -- convex set

### Output

If the matrix is null, a `ZeroSet` is returned; otherwise a lazy linear map.
"""
*(M::AbstractMatrix, X::LazySet) = LinearMap(M, X)


"""
```
    *(a::N, X::LazySet) where {N}
```

Return a linear map of a convex set by a scalar value.

### Input

- `a` -- scalar
- `X` -- convex set

### Output

The linear map of the convex set.
"""
function *(a::N, X::LazySet) where {N}
    return LinearMap(a * speye(N, dim(X)), X)
end

"""
```
    *(a::N, map::S)::S where {S<:LinearMap{NM, N}} where {NM, N}
```

Return a linear map scaled by a scalar value.

### Input

- `a`   -- scalar
- `map` -- linear map

### Output

The scaled linear map.
"""
function *(a::N, map::S)::S where {S<:LinearMap{NM, N}} where {NM, N}
    return LinearMap(a * map.M, map.X)
end

"""
```
    *(M::AbstractMatrix, Z::ZeroSet)
````
Linear map of a zero set.

### Input

- `M` -- abstract matrix
- `Z` -- zero set

### Output

The zero set with the output dimension of the linear map.
"""
function *(M::AbstractMatrix, Z::ZeroSet)
    @assert dim(Z) == size(M, 2)
    return ZeroSet(size(M, 1))
end

"""
    dim(lm::LinearMap)::Int

Return the dimension of a linear map.

### Input

- `lm` -- linear map

### Output

The ambient dimension of the linear map.
"""
function dim(lm::LinearMap)::Int
    return size(lm.M, 1)
end

"""
    σ(d::AbstractVector{<:Real}, lm::LinearMap)::AbstractVector{<:Real}

Return the support vector of the linear map.

### Input

- `d`  -- direction
- `lm` -- linear map

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``L = M⋅S``, where ``M`` is a matrix and ``S`` is a convex set, it follows
that ``σ(d, L) = M⋅σ(M^T d, S)`` for any direction ``d``.
"""
function σ(d::V, lm::LinearMap) where {N<:Real, V<:AbstractVector{N}}
    return lm.M * σ(lm.M.' * d, lm.X)
end

"""
    ∈(x::AbstractVector{N}, lm::LinearMap{NM, N})::Bool where {NM, N<:Real}

Check whether a given point is contained in a linear map of a convex set.

### Input

- `x`  -- point/vector
- `lm` -- linear map of a convex set

### Output

`true` iff ``x ∈ lm``.

### Algorithm

Note that ``x ∈ M⋅S`` iff ``M^{-1}⋅x ∈ S``.
This implementation does not explicitly invert the matrix, which is why it also
works for non-square matrices.

### Examples

```jldoctest
julia> lm = LinearMap([2.0 0.0; 0.0 1.0], BallInf([1., 1.], 1.));

julia> ∈([5.0, 1.0], lm)
false
julia> ∈([3.0, 1.0], lm)
true
```

An example with non-square matrix:
```jldoctest
julia> B = BallInf(zeros(4), 1.);

julia> M = [1. 0 0 0; 0 1 0 0]/2;

julia> ∈([0.5, 0.5], M*B)
true
```
"""
function ∈(x::AbstractVector{N}, lm::LinearMap{NM, N})::Bool where {NM, N<:Real}
    return ∈(lm.M \ x, lm.X)
end
