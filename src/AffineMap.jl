import Base: isempty

export AffineMap,
       an_element,
       isempty

"""
    AffineMap{N<:Real, S<:LazySet{N},
             NM, MAT<:AbstractMatrix{NM}, VN<:AbstractVector{NM}} <: LazySet{N}

Type that represents an affine transformation ``M⋅X ⊕ b`` of a convex set ``X``.

### Fields

- `M` -- matrix/linear map
- `b` -- vector
- `X` -- convex set

### Notes

The affine map is the composition of a linear map and a translation.
This type is parametric in the coefficients of the linear map, `NM`, which may be
different from the numeric type of the wrapped set (`N`).
"""
struct AffineMap{N<:Real, S<:LazySet{N},
                 NM, MAT<:AbstractMatrix{NM}, VN<:AbstractVector{NM}} <: LazySet{N}
    M::MAT
    b::VN
    X::S

    # default constructor with dimension match check
    function AffineMap{N, S, NM, MAT, VN}(M::MAT, b::VN, X::S) where {N<:Real,
                        S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM}, VN}

        @assert dim(X) == size(M, 2) "a matrix of size $(size(M)) cannot be "
            "applied to a set of dimension $(dim(X))"

        @assert size(M, 1) == length(b) "the output dimension of the map, $(size(M, 2)), "
            "is incompatible with the dimension of the translation vector, $(dim(b))"

        return new{N, S, NM, MAT, VN}(M, b, X)
    end
end

# ============================ 
# LazySet interface functions
# ============================

"""
    dim(am::AffineMap)::Int

Return the dimension of an affine map.

### Input

- `am` -- affine map

### Output

The dimension of an affine map.
"""
function dim(am::AffineMap)::Int
    return length(am.b)
end


"""
    σ(d::AbstractVector{N}, am::AffineMap{N}) where {N<:Real}

Return the support vector of an affine map.

### Input

- `d`  -- direction
- `am` -- affine map

### Output

The support vector in the given direction.
"""
function σ(d::AbstractVector{N}, am::AffineMap{N}) where {N<:Real}
    return am.M * σ(_At_mul_B(am.M, d), am.X) + am.b
end

"""
    ρ(d::AbstractVector{N}, am::AffineMap{N}) where {N<:Real}

Return the support function of an affine map.

### Input

- `d`  -- direction
- `am` -- affine map

### Output

The support function in the given direction.
"""
function ρ(d::AbstractVector{N}, am::AffineMap{N}) where {N<:Real}
    return ρ(_At_mul_B(am.M, d), am.X) + dot(d, am.b)
end

"""
    an_element(am::AffineMap)

Return some element of an affine map.

### Input

- `am` -- affine map

### Output

An element of the affine map. It relies on the `an_element` function of the
wrapped set.
"""
function an_element(am::AffineMap)
    return am.M * an_element(am.X) + am.b
end

"""
    isempty(am::AffineMap)::Bool

Return whether an affine map is empty or not.

### Input

- `am` -- affine map

### Output

`true` iff the wrapped set is empty and the affine vector is empty.
"""
function isempty(am::AffineMap)::Bool
    return isempty(am.X) && isempty(am.b)
end
