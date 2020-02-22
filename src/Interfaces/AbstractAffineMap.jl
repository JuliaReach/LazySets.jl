import Base: isempty, ∈

export AbstractAffineMap,
       matrix, vector, set

"""
    AbstractAffineMap{N<:Real, S<:LazySet{N}} <: LazySet{N}

Abstract type for affine maps.

### Notes

See [`AffineMap`](@ref) for a standard implementation of this interface.

Every concrete `AbstractAffineMap` must define the following functions:
- `matrix(::AbstractAffineMap)` -- return the linear map
- `vector(::AbstractAffineMap)` -- return the affine translation vector
- `set(::AbstractAffineMap)` -- return the set that the map is applied to

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractAffineMap)
6-element Array{Any,1}:
 AffineMap
 ExponentialMap
 ExponentialProjectionMap
 LinearMap
 ResetMap
 Translation
```
"""
abstract type AbstractAffineMap{N<:Real, S<:LazySet{N}} <: LazySet{N} end

isoperationtype(::Type{<:AbstractAffineMap}) = true
isconvextype(::Type{<:AbstractAffineMap{N, S}}) where {N, S} = isconvextype(S)


# --- common AbstractAffineMap functions ---


"""
    dim(am::AbstractAffineMap)

Return the dimension of an affine map.

### Input

- `am` -- affine map

### Output

The dimension of an affine map.
"""
function dim(am::AbstractAffineMap)
    return length(vector(am))
end

"""
    σ(d::AbstractVector{N}, am::AbstractAffineMap{N}) where {N<:Real}

Return the support vector of an affine map.

### Input

- `d`  -- direction
- `am` -- affine map

### Output

The support vector in the given direction.
"""
function σ(d::AbstractVector{N}, am::AbstractAffineMap{N}) where {N<:Real}
    A = matrix(am)
    return A * σ(_At_mul_B(A, d), set(am)) + vector(am)
end

"""
    ρ(d::AbstractVector{N}, am::AbstractAffineMap{N}) where {N<:Real}

Return the support function of an affine map.

### Input

- `d`  -- direction
- `am` -- affine map

### Output

The support function in the given direction.
"""
function ρ(d::AbstractVector{N}, am::AbstractAffineMap{N}) where {N<:Real}
    return ρ(_At_mul_B(matrix(am), d), set(am)) + dot(d, vector(am))
end

"""
    an_element(am::AbstractAffineMap)

Return some element of an affine map.

### Input

- `am` -- affine map

### Output

An element of the affine map. It relies on the `an_element` function of the
wrapped set.
"""
function an_element(am::AbstractAffineMap)
    return matrix(am) * an_element(set(am)) + vector(am)
end

"""
    isempty(am::AbstractAffineMap)

Return whether an affine map is empty or not.

### Input

- `am` -- affine map

### Output

`true` iff the wrapped set is empty.
"""
function isempty(am::AbstractAffineMap)
    return isempty(set(am))
end

"""
    isbounded(am::AbstractAffineMap; cond_tol::Number=DEFAULT_COND_TOL)

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
function isbounded(am::AbstractAffineMap; cond_tol::Number=DEFAULT_COND_TOL)
    M = matrix(am)
    if iszero(M) || isbounded(set(am))
        return true
    end
    if isinvertible(M; cond_tol=cond_tol)
        return false
    end
    return isbounded_unit_dimensions(am)
end

"""
    ∈(x::AbstractVector{N}, am::AbstractAffineMap{N}) where {N<:Real}

Check whether a given point is contained in the affine map of a convex set.

### Input

- `x`  -- point/vector
- `am` -- affine map of a convex set

### Output

`true` iff ``x ∈ am``.

### Algorithm

Note that ``x ∈ M⋅S ⊕ v`` iff ``M^{-1}⋅(x - v) ∈ S``.
This implementation does not explicitly invert the matrix, which is why it also
works for non-square matrices.

### Examples

```jldoctest
julia> am = AffineMap([2.0 0.0; 0.0 1.0], BallInf([1., 1.], 1.), [-1.0, -1.0]);

julia> [5.0, 1.0] ∈ am
false

julia> [3.0, 1.0] ∈ am
true
```

An example with a non-square matrix:

```jldoctest
julia> B = BallInf(zeros(4), 1.);

julia> M = [1. 0 0 0; 0 1 0 0]/2;

julia> [0.5, 0.5] ∈ M*B
true
```
"""
function ∈(x::AbstractVector{N}, am::AbstractAffineMap{N}) where {N<:Real}
    return matrix(am) \ (x - vector(am)) ∈ set(am)
end

"""
    vertices_list(am::AbstractAffineMap{N};
                  [apply_convex_hull]::Bool) where {N<:Real}

Return the list of vertices of a (polyhedral) affine map.

### Input

- `am`                -- affine map
- `apply_convex_hull` -- (optional, default: `true`) if `true`, apply the convex
                         hull operation to the list of vertices transformed by
                         the affine map

### Output

A list of vertices.

### Algorithm

This implementation computes all vertices of `X`, then transforms them through
the affine map, i.e. `x ↦ M*x + v` for each vertex `x` of `X`. By default, the
convex hull operation is taken before returning this list. For dimensions three
or higher, this operation relies on the functionality through the concrete
polyhedra library `Polyhedra.jl`.

If you are not interested in taking the convex hull of the resulting vertices
under the affine map, pass `apply_convex_hull=false` as a keyword argument.

Note that we assume that the underlying set `X` is polyhedral, either concretely
or lazily, i.e. there the function `vertices_list` should be applicable.
"""
function vertices_list(am::AbstractAffineMap{N};
                       apply_convex_hull::Bool=true) where {N<:Real}
    # for a zero linear map, the result is just the affine translation
    A = matrix(am)
    b = vector(am)
    if iszero(A)
        return [b]
    end

    # collect vertices list of the wrapped set
    X = set(am)
    @assert applicable(vertices_list, X)  "this function requires that the " *
        "list of vertices is available, but it is not"
    vlist_X = vertices_list(X)

    # create resulting vertices list
    vlist = [A * x + b for x in vlist_X]

    return apply_convex_hull ? convex_hull!(vlist) : vlist
end

"""
    constraints_list(am::AbstractAffineMap{N}) where {N<:Real}

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
function constraints_list(am::AbstractAffineMap{N}) where {N<:Real}
    return _constraints_list_translation(LinearMap(matrix(am), set(am)),
                                         vector(am))
end

"""
    linear_map(M::AbstractMatrix{N}, am::AbstractAffineMap{N}) where {N<:Real}

Return the linear map of a lazy affine map.

### Input

- `M`  -- matrix
- `am` -- affine map

### Output

A set corresponding to the linear map of the lazy affine map of a set.
"""
function linear_map(M::AbstractMatrix{N},
                    am::AbstractAffineMap{N}) where {N<:Real}
     return translate(linear_map(M * matrix(am), set(am)), M * vector(am))
end
