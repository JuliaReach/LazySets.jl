export AbstractAffineMap,
       matrix, vector, set

"""
    AbstractAffineMap{N, S<:LazySet{N}} <: LazySet{N}

Abstract type for affine maps.

### Notes

See [`AffineMap`](@ref) for a standard implementation of this interface.

Every concrete `AbstractAffineMap` must define the following methods:

- `matrix(::AbstractAffineMap)` -- return the linear map
- `vector(::AbstractAffineMap)` -- return the affine translation vector
- `set(::AbstractAffineMap)` -- return the set that the map is applied to

The subtypes of `AbstractAffineMap`:

```jldoctest; setup = :(using LazySets: subtypes)
julia> subtypes(AbstractAffineMap)
7-element Vector{Any}:
 AffineMap
 ExponentialMap
 ExponentialProjectionMap
 InverseLinearMap
 LinearMap
 ResetMap
 Translation
```
"""
abstract type AbstractAffineMap{N,S<:LazySet{N}} <: LazySet{N} end

"""
    matrix(X::AbstractAffineMap)

Return the matrix of an affine map.

### Input

- `X` -- affine map

### Output

The matrix of `X`.
"""
function matrix(::AbstractAffineMap) end

"""
    vector(X::AbstractAffineMap)

Return the vector of an affine map.

### Input

- `X` -- affine map

### Output

The vector of `X`.
"""
function vector(::AbstractAffineMap) end

"""
    set(X::AbstractAffineMap)

Return the set of an affine map.

### Input

- `X` -- affine map

### Output

The set of `X` before applying the map.
"""
function set(::AbstractAffineMap) end

isoperationtype(::Type{<:AbstractAffineMap}) = true

isconvextype(::Type{<:AbstractAffineMap{N,S}}) where {N,S} = isconvextype(S)

ispolyhedraltype(::Type{<:AbstractAffineMap{N,S}}) where {N,S} = ispolyhedraltype(S)

function ispolyhedral(am::AbstractAffineMap)
    ispolyhedraltype(typeof(am)) && return true

    return ispolyhedral(set(am))
end

function dim(am::AbstractAffineMap)
    return length(vector(am))
end

@validate function σ(d::AbstractVector, am::AbstractAffineMap)
    A = matrix(am)
    return A * σ(At_mul_B(A, d), set(am)) + vector(am)
end

@validate function ρ(d::AbstractVector, am::AbstractAffineMap)
    return ρ(At_mul_B(matrix(am), d), set(am)) + dot(d, vector(am))
end

"""
# Extended help

### Algorithm

The implementation relies on the `an_element` method of the wrapped set.
"""
function an_element(am::AbstractAffineMap)
    return matrix(am) * an_element(set(am)) + vector(am)
end

"""
# Extended help

    isempty(am::AbstractAffineMap)

### Algorithm

The result is `true` iff the wrapped set is empty.
"""
function isempty(am::AbstractAffineMap)
    return isempty(set(am))
end

"""
# Extended help

    isbounded(am::AbstractAffineMap; cond_tol::Number=DEFAULT_COND_TOL)

### Input

- `cond_tol` -- (optional) tolerance of matrix condition (used to check whether
                the matrix is invertible)

### Algorithm

We first check if the matrix is zero or the wrapped set is bounded.
If not, we perform a sufficient check whether the matrix is invertible.
If the matrix is invertible, then the map being bounded is equivalent to the
wrapped set being bounded, and hence the map is unbounded.
Otherwise, we check boundedness via [`_isbounded_unit_dimensions`](@ref).
"""
function isbounded(am::AbstractAffineMap; cond_tol::Number=DEFAULT_COND_TOL)
    M = matrix(am)
    if iszero(M) || isbounded(set(am))
        return true
    end
    if isinvertible(M; cond_tol=cond_tol)
        return false
    end
    return _isbounded_unit_dimensions(am)
end

"""
# Extended help

    in(x::AbstractVector, am::AbstractAffineMap)

### Algorithm

Observe that ``x ∈ M⋅S ⊕ v`` iff ``M^{-1}⋅(x - v) ∈ S``.
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
@validate function in(x::AbstractVector, am::AbstractAffineMap)
    if !iswellconditioned(matrix(am))
        # ill-conditioned matrix; use concrete set representation
        return x ∈ affine_map(matrix(am), set(am), vector(am))
    end
    return matrix(am) \ (x - vector(am)) ∈ set(am)
end

"""
# Extended help

    center(am::AbstractAffineMap)

### Algorithm

The implementation relies on the `center` method of the wrapped set.
"""
function center(am::AbstractAffineMap)
    return matrix(am) * center(set(am)) + vector(am)
end

"""
# Extended help

    vertices_list(am::AbstractAffineMap; [apply_convex_hull]::Bool=true)

### Input

- `apply_convex_hull` -- (optional, default: `true`) if `true`, apply the convex
                         hull operation to the list of vertices transformed by
                         the affine map

### Algorithm

This implementation computes all vertices of `X`, then transforms them through
the affine map, i.e., `x ↦ M*x + v` for each vertex `x` of `X`. By default, the
convex-hull operation is taken before returning this list. For dimensions three
or higher, this operation relies on the functionality through the concrete
polyhedra library `Polyhedra.jl`.

If you are not interested in taking the convex hull of the resulting vertices
under the affine map, pass `apply_convex_hull=false` as a keyword argument.

Note that we assume that the underlying set `X` is polytopic, either concretely
or lazily, i.e., the function `vertices_list` should be applicable.
"""
@validate function vertices_list(am::AbstractAffineMap; apply_convex_hull::Bool=true)
    # for a zero linear map, the result is just the affine translation
    A = matrix(am)
    b = vector(am)
    if iszero(A)
        return [b]
    end

    # collect vertices list of the wrapped set
    X = set(am)
    vlist_X = vertices_list(X)

    # create resulting vertices list
    vlist = [A * x + b for x in vlist_X]

    return apply_convex_hull ? convex_hull!(vlist) : vlist
end

"""
# Extended help

    constraints_list(am::AbstractAffineMap)

### Notes

We assume that the underlying set `X` is polyhedral, i.e., offers a method
`constraints_list(X)`.

### Algorithm

This implementation uses the `constraints_list` method to compute the list of
constraints of the translation of a lazy [`LinearMap`](@ref).
"""
function constraints_list(am::AbstractAffineMap)
    return _constraints_list_translation(LinearMap(matrix(am), set(am)),
                                         vector(am))
end

"""
# Extended help

    isuniversal(am::AbstractAffineMap)

### Algorithm

An affine map is universal iff the wrapped set is universal and the matrix does
not map any dimension to zero.
"""
function isuniversal(am::AbstractAffineMap)
    if isuniversal(set(am))
        for row in eachrow(matrix(am))
            if all(isapproxzero, row)
                return false
            end
        end
        return true
    end
    return false
end

@validate function linear_map(M::AbstractMatrix, am::AbstractAffineMap)
    return affine_map(M * matrix(am), set(am), M * vector(am))
end

"""
# Extended help

    volume(am::AbstractAffineMap)

### Notes

This implementation requires a dimension-preserving map (i.e., a square matrix).

### Algorithm

A square linear map scales the volume of any set by its absolute determinant.
A translation does not affect the volume.
Thus, the volume of `M * X + {v}` is `|det(M)| * volume(X)`.
"""
function volume(am::AbstractAffineMap)
    checksquare(matrix(am))

    return abs(det(matrix(am))) * volume(set(am))
end
