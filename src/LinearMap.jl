import Base: *, ∈

export LinearMap,
       an_element

"""
    LinearMap{N<:Real, S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM}} <: LazySet{N}

Type that represents a linear transformation ``M⋅S`` of a convex set ``S``.

### Fields

- `M` -- matrix/linear map
- `X` -- convex set

### Notes

This type is parametric in the elements of the linear map, `NM`, which is
independent of the numeric type of the target set (`N`).
Typically `NM = N`, but there may be exceptions, e.g., if `NM` is an interval
that holds numbers of type `N`, where `N` is a floating point number type such
as `Float64`.
"""
struct LinearMap{N<:Real, S<:LazySet{N},
                 NM, MAT<:AbstractMatrix{NM}} <: LazySet{N}
    M::MAT
    X::S

    # default constructor with dimension match check
    function LinearMap{N, S, NM, MAT}(M::MAT, X::S) where {N<:Real,
            S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM}}
        @assert dim(X) == size(M, 2) "a linear map of size $(size(M)) cannot " *
            "be applied to a set of dimension $(dim(X))"
        return new{N, S, NM, MAT}(M, X)
    end
end

# convenience constructor without type parameter
LinearMap(M::MAT,
          X::S) where {N<:Real, S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM}} =
    LinearMap{N, S, NM, MAT}(M, X)

# constructor from a linear map: perform the matrix multiplication immediately
LinearMap(M::MAT, lm::LinearMap{N, S, NM, MAT}) where {N<:Real, S<:LazySet{N},
        NM, MAT<:AbstractMatrix{NM}} =
    LinearMap{N, S, NM, MAT}(M * lm.M, lm.X)

"""
```
    *(M::AbstractMatrix, X::LazySet)
```

Return the linear map of a convex set.

### Input

- `M` -- matrix/linear map
- `X` -- convex set

### Output

A lazy linear map, i.e. a `LinearMap` instance.
"""
*(M::AbstractMatrix, X::LazySet) = LinearMap(M, X)

*(M::AbstractVector, X::LazySet) = LinearMap(reshape(M, length(M), 1), X)

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
function *(a::N, X::LazySet{N}) where {N}
    n = dim(X)
    return LinearMap(sparse(a * I, n, n), X)
end

"""
```
    *(a::N, lm::LM)::LM where {N<:Real, S<:LazySet{N}, LM<:LinearMap{N, S}}
```

Return a linear map scaled by a scalar value.

### Input

- `a`  -- scalar
- `lm` -- linear map

### Output

The scaled linear map.
"""
function *(a::N, lm::LM)::LM where {N<:Real, S<:LazySet{N}, LM<:LinearMap{N, S}}
    return LinearMap(a * lm.M, lm.X)
end

"""
```
    *(M::MAT, Z::ZeroSet{N})::ZeroSet{N} where {N<:Real, MAT<:AbstractMatrix}
```

A linear map of a zero set, which is simplified to a zero set (the absorbing
element).

### Input

- `M` -- abstract matrix
- `Z` -- zero set

### Output

The zero set with the output dimension of the linear map.
"""
function *(M::MAT,
           Z::ZeroSet{N})::ZeroSet{N} where {N<:Real, MAT<:AbstractMatrix}
    @assert dim(Z) == size(M, 2) "a linear map of size $(size(M)) cannot " *
            "be applied to a set of dimension $(dim(Z))"
    return ZeroSet{N}(size(M, 1))
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
    σ(d::AbstractVector{N}, lm::LinearMap{N}) where {N<:Real}

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
function σ(d::AbstractVector{N}, lm::LinearMap{N}) where {N<:Real}
    return lm.M * σ(_At_mul_B(lm.M, d), lm.X)
end

"""
    ρ(d::AbstractVector{N}, lm::LinearMap{N}) where {N<:Real}

Return the support function of the linear map.

### Input

- `d`  -- direction
- `lm` -- linear map

### Output

The support function in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``L = M⋅S``, where ``M`` is a matrix and ``S`` is a convex set, it follows
that ``ρ(d, L) = ρ(M^T d, S)`` for any direction ``d``.
"""
function ρ(d::AbstractVector{N}, lm::LinearMap{N}) where {N<:Real}
    return ρ(_At_mul_B(lm.M, d), lm.X)
end

"""
    ∈(x::AbstractVector{N}, lm::LinearMap{N})::Bool where {N<:Real}

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
function ∈(x::AbstractVector{N}, lm::LinearMap{N})::Bool where {N<:Real}
    return ∈(lm.M \ x, lm.X)
end

"""
    an_element(lm::LinearMap)

Return some element of a linear map.

### Input

- `lm` -- linear map

### Output

An element in the linear map.
It relies on the `an_element` function of the wrapped set.
"""
function an_element(lm::LinearMap)
    return lm.M * an_element(lm.X)
end

"""
    vertices_list(lm::LinearMap{N})::Vector{Vector{N}} where N<:Real

Return the list of vertices of a (polytopic) linear map.

### Input

- `lm` -- linear map

### Output

A list of vertices.

### Algorithm

We assume that the underlying set `X` is polytopic.
Then the result is just the linear map applied to the vertices of `X`.
"""
function vertices_list(lm::LinearMap{N})::Vector{Vector{N}} where N<:Real
    # for a zero map, the result is just the list containing the origin
    if iszero(lm.M)
        return [zeros(N, dim(lm))]
    end

    # collect low-dimensional vertices lists
    vlist_X = vertices_list(lm.X)

    # create resulting vertices list
    vlist = Vector{Vector{N}}()
    sizehint!(vlist, length(vlist_X))
    for v in vlist_X
        push!(vlist, lm.M * v)
    end

    return vlist
end
