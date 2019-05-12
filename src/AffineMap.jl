import Base: isempty

export AffineMap,
       an_element,
       isempty,
       isbounded,
       ∈,
       vertices_list,
       constraints_list,
       linear_map

"""
    AffineMap{N<:Real, S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM},
             VN<:AbstractVector{NM}} <: LazySet{N}

Type that represents an affine transformation ``M⋅X ⊕ b`` of a convex set ``X``.

### Fields

- `M` -- matrix/linear map
- `b` -- vector
- `X` -- convex set

### Notes

The affine map is the composition of a linear map and a translation. This type is
parametric in the coefficients of the linear map, `NM`, which may be different from
the numeric type of the wrapped set (`N`). However, the numeric type of the
translation vector should be `NM`.
"""
struct AffineMap{N<:Real, S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM},
                VN<:AbstractVector{NM}} <: LazySet{N}
    M::MAT
    b::VN
    X::S

    # default constructor with dimension match check
    function AffineMap{N, S, NM, MAT, VN}(M::MAT, b::VN, X::S) where {N<:Real,
                        S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM}, VN}

        @assert dim(X) == size(M, 2) "a matrix of size $(size(M)) cannot be " *
            "applied to a set of dimension $(dim(X))"

        @assert size(M, 1) == length(b) "a map with output dimension $(size(M, 1)) " *
            "is incompatible with the dimension of the translation vector, $(length(b))"

        return new{N, S, NM, MAT, VN}(M, b, X)
    end
end

# convenience constructor without type parameter
AffineMap(M::MAT, b::VN, X::S) where {N<:Real, S<:LazySet{N}, NM,
    MAT<:AbstractMatrix{NM}, VN<:AbstractVector{NM}} = AffineMap{N, S, NM, MAT, VN}(M, b, X)

# ============================ 
# Arithmetic functions
# ============================

"""
```
    *(M::AbstractMatrix{N}, am::AffineMap{N}) where {N<:Real}
```

Transform an affine map under matrix multiplication.

### Input

- `M`  -- matrix
- `am` -- affine map

### Output

A lazy affine map, i.e. an `AffineMap`, such that the new map is related to the
old one through `am.M ↦ M * am.M` and `am.b ↦ M * am.b`.
"""
function *(M::AbstractMatrix{N}, am::AffineMap{N}) where {N<:Real}
    return AffineMap(M * am.M, M * am.b, am.X)
end

"""
```
    *(α::N, am::AffineMap{N}) where {N<:Real}
```

Return the affine map scaled by a given number.

### Input

- `α`  -- scalar
- `am` -- affine map

### Output

The scaled affine map.
"""
function *(α::N, am::AffineMap{N}) where {N<:Real}
    return AffineMap(α * am.M, α * am.b, am.X)
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

"""
    isbounded(am::AffineMap; cond_tol::Number=DEFAULT_COND_TOL)::Bool

Determine whether an affine map is bounded.

### Input

- `am`       -- affine map
- `cond_tol` -- (optional) tolerance of matrix condition (used to check whether
                the matrix is invertible)

### Output

`true` iff the affine map is bounded.

### Algorithm

We first check if the matrix is zero or the wrapped set is bounded.
If not, we perform a sufficient check whether the matrix is invertible.
If the matrix is invertible, then the map being bounded is equivalent to the
wrapped set being bounded, and hence the map is unbounded.
Otherwise, we check boundedness via [`isbounded_unit_dimensions`](@ref).
"""
function isbounded(am::AffineMap; cond_tol::Number=DEFAULT_COND_TOL)::Bool
    if iszero(am.M) || isbounded(am.X)
        return true
    end
    if isinvertible(am.M; cond_tol=cond_tol)
        return false
    end
    return isbounded_unit_dimensions(am)
end

"""
    ∈(x::AbstractVector{N}, am::AffineMap{N})::Bool where {N<:Real}

Check whether a given point is contained in the affine map of a convex set.

### Input

- `x`  -- point/vector
- `am` -- affine map of a convex set

### Output

`true` iff ``x ∈ am``.

### Algorithm

Note that ``x ∈ M⋅S ⊕ b`` iff ``M^{-1}⋅(x - b) ∈ S``.
This implementation does not explicitly invert the matrix, which is why it also
works for non-square matrices.

### Examples

```jldoctest
julia> am = AffineMap([2.0 0.0; 0.0 1.0], [-1.0, -1.0], BallInf([1., 1.], 1.));

julia> ∈([5.0, 1.0], am)
false

julia> ∈([3.0, 1.0], am)
true
```

An example with a non-square matrix:

```jldoctest
julia> B = BallInf(zeros(4), 1.);

julia> M = [1. 0 0 0; 0 1 0 0]/2;

julia> ∈([0.5, 0.5], M*B)
true
```
"""
function ∈(x::AbstractVector{N}, am::AffineMap{N})::Bool where {N<:Real}
    return ∈(am.M \ (x - am.b), am.X)
end

"""
    vertices_list(am::AffineMap{N};
                  [apply_convex_hull]::Bool)::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a (polyhedral) affine map.

### Input

- `am`                -- affine map
- `apply_convex_hull` -- (optional, default: `true`) if `true`, apply the convex
                         hull operation to the list of vertices transformed by the
                         affine map 

### Output

A list of vertices.

### Algorithm

This implementation computes all vertices of `X`, then transforms them through
the affine map, i.e.  `v ↦ M*v + b` for each vertex `v` of `X`. By default, the
convex hull operation is taken before returning this list. For dimensions three
or higher, this operation relies on the functionality through the concrete
polyhedra library `Polyhedra.jl`.

If you are not interested in taking the convex hull of the resulting vertices under
the affine map, pass `apply_convex_hull=false` as a keyword argument.

Note that we assume that the underlying set `X` is polyhedral, either concretely
or lazily, i.e. there the function `vertices_list` should be applicable.
"""
function vertices_list(am::AffineMap{N};
                       apply_convex_hull::Bool=true)::Vector{Vector{N}} where {N<:Real}
    # for a zero map, the result is just the singleton containing the translation
    if iszero(am.M)
        return Singleton(am.b)
    end

    # collect low-dimensional vertices lists
    vlist_X = vertices_list(am.X)

    # create resulting vertices list
    vlist = Vector{Vector{N}}()
    sizehint!(vlist, length(vlist_X))
    for v in vlist_X
        push!(vlist, am.M * v + am.b)
    end

    return apply_convex_hull ? convex_hull(vlist) : vlist
end

"""
    constraints_list(am::AffineMap{N}) where {N<:Real}

Return the list of constraints of a (polyhedral) affine map.

### Input

- `am` -- affine map

### Output

The list of constraints of the affine map.

### Notes

We assume that the underlying set `X` is polyhedral, i.e., offers a method
`constraints_list(X)`.

### Algorithm

Falls back to the list of constraints of the translation of a lazy linear map.
"""
function constraints_list(am::AffineMap{N}) where {N<:Real}
    return constraints_list(LinearMap(am.M, am.X) ⊕ am.b)
end

"""
    linear_map(M::AbstractMatrix{N}, am::AffineMap{N}) where {N<:Real}

Return the linear map of a lazy affine map.

### Input

- `M`  -- matrix
- `am` -- affine map

### Output

The polytope representing the linear map of the lazy affine map of a set.  
"""
function linear_map(M::AbstractMatrix{N}, am::AffineMap{N}) where {N<:Real}
     return linear_map(M * am.M, am.X) + M * am.b
end
