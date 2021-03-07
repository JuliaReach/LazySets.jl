import Base: *, ∈, isempty

export InverseLinearMap,
       an_element,
       constraints_list

"""
   InverseLinearMap{N, S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM}} <: AbstractAffineMap{N, S}

Given a linear transformation ``M``, this type represents the linear
transformation ``M⁻¹⋅X`` of a set ``X`` without actually computing ``M⁻¹``.

### Fields

- `M` -- invertible matrix
- `X` -- set

### Notes

This type is parametric in the elements of the inverse linear map, `NM`, which is
independent of the numeric type of the wrapped set (`N`). Typically `NM = N`,
but there may be exceptions, e.g., if `NM` is an interval that holds numbers
of type `N`, where `N` is a floating point number type such as `Float64`.

### Examples

For the examples we create a ``3×3`` matrix and a unit three-dimensional square.

```jldoctest ilp_constructor
julia> A = [1 2 3; 2 3 1; 3 1 2];

julia> X = BallInf([0, 0, 0], 1);

julia> ilm = InverseLinearMap(A, X)
InverseLinearMap{Int64,BallInf{Int64,Array{Int64,1}},Int64,Array{Int64,2}}([1 2 3; 2 3 1; 3 1 2], BallInf{Int64,Array{Int64,1}}([0, 0, 0], 1))
```

Applying an inverse linear map to a `InverseLinearMap` object combines the two maps into
a single `InverseLinearMap` instance.

```jldoctest ilp_constructor
julia> B = transpose(A); ilm2 = InverseLinearMap(B, ilm)
InverseLinearMap{Int64,BallInf{Int64,Array{Int64,1}},Int64,Array{Int64,2}}([14 11 11; 11 14 11; 11 11 14], BallInf{Int64,Array{Int64,1}}([0, 0, 0], 1))

julia> ilm2.M == A*B
true
```

The application of an `InverseLinearMap` to a `ZeroSet` or an `EmptySet` is
simplified automatically.

```jldoctest ilp_constructor
julia> InverseLinearMap(A, ZeroSet{Int}(3))
ZeroSet{Int64}(3)
```
"""
struct InverseLinearMap{N, S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM}} <: AbstractAffineMap{N, S}
    M::MAT
    X::S

    # default constructor with dimension match check
    function InverseLinearMap(M::MAT, X::S; check_invertibility::Bool=false) where {N, S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM}}
        @assert dim(X) == size(M, 1) "a linear map of size $(size(M)) cannot " *
                                     "be applied to a set of dimension $(dim(X))"
        if check_invertibility
            @assert isinvertible(M) "the linear map is not invertible"
        end
        return new{N, S, NM, MAT}(M, X)
    end
end

# convenience constructor from a UniformScaling
function InverseLinearMap(M::UniformScaling{N}, X::LazySet) where {N}
    return InverseLinearMap(M.λ, X)
end

# convenience constructor from a scalar
function InverseLinearMap(α::N, X::LazySet) where {N<:Real}
    if α == one(N)
        return X
    end
    D = Diagonal(fill(α, dim(X)))
    return InverseLinearMap(D, X)
end

# combine two linear maps into a single linear map
function InverseLinearMap(M::AbstractMatrix, ilm::InverseLinearMap)
    return InverseLinearMap(ilm.M * M, ilm.X)
end

# ZeroSet is "almost absorbing" for InverseLinearMap (only the dimension changes)
function InverseLinearMap(M::AbstractMatrix{N}, Z::ZeroSet{N}) where {N}
    @assert dim(Z) == size(M, 2) "a linear map of size $(size(M)) cannot " *
            "be applied to a set of dimension $(dim(Z))"
    return ZeroSet{N}(size(M, 1))
end

# EmptySet is absorbing for LinearMap
function InverseLinearMap(M::AbstractMatrix, ∅::EmptySet)
    return ∅
end

# --- AbstractAffineMap interface functions ---

function matrix(ilm::InverseLinearMap)
    return inv(ilm.M)
end

function vector(ilm::InverseLinearMap{N}) where {N}
    return spzeros(N, dim(ilm))
end

function set(ilm::InverseLinearMap)
    return ilm.X
end

"""
    dim(ilm::InverseLinearMap)

Return the dimension of an inverse linear map.

### Input

- `ilm` -- inverse linear map

### Output

The ambient dimension of the inverse linear map.
"""
function dim(ilm::InverseLinearMap)
    return size(ilm.M, 1)
end

"""
    σ(d::AbstractVector, ilm::InverseLinearMap)

Return the support vector of the inverse linear map.

### Input

- `d`   -- direction
- `ilm` -- inverse linear map

### Output

The support vector in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``L = M^{-1}⋅X``, where ``M`` is a matrix and ``X`` is a set, since
(M^T)^{-1}=(M^{-1})^T, it follows that ``σ(d, L) = M^{-1}⋅σ((M^T)^{-1} d, X)``
for any direction ``d``.
"""
function σ(d::AbstractVector, ilm::InverseLinearMap)
    y = transpose(ilm.M) \ d
    return ilm.M \ σ(y, ilm.X)
end

"""
    ρ(d::AbstractVector, ilm::InverseLinearMap)

Return the support function of the inverse linear map.

### Input

- `d`      -- direction
- `ilm`    -- inverse linear map

### Output

The support function in the given direction.
If the direction has norm zero, the result depends on the wrapped set.

### Notes

If ``L = M^{-1}⋅X``, where ``M`` is a matrix and ``X`` is a set, it follows
that ``ρ(d, L) = ρ((M^T)^{-1} d, X)`` for any direction ``d``.
"""
function ρ(d::AbstractVector, ilm::InverseLinearMap)
    y = transpose(ilm.M) \ d
    return ρ(y, ilm.X)
end

"""
    ∈(x::AbstractVector, ilm::InverseLinearMap)

Check whether a given point is contained in the inverse linear map of a set.

### Input

- `x`   -- point/vector
- `ilm` -- inverse linear map of a set

### Output

`true` iff ``x ∈ ilm``.

### Algorithm

This implementation does not explicitly invert the matrix since it uses the
property ``x ∈ M^{-1}⋅X`` iff ``M⋅x ∈ X``.

### Examples

```jldoctest
julia> ilm = LinearMap([0.5 0.0; 0.0 -0.5], BallInf([0., 0.], 1.));

julia> [1.0, 1.0] ∈ ilm
false

julia> [0.1, 0.1] ∈ ilm
true
```
"""
function ∈(x::AbstractVector, ilm::InverseLinearMap)
    y = ilm.M * x
    return y ∈ ilm.X
end

"""
    an_element(ilm::InverseLinearMap)

Return some element of an inverse linear map.

### Input

- `ilm` -- inverse linear map

### Output

An element in the inverse linear map.
It relies on the `an_element` function of the wrapped set.
"""
function an_element(lm::InverseLinearMap)
    return lm.M \ an_element(lm.X)
end

"""
    vertices_list(ilm::InverseLinearMap; prune::Bool=true)

Return the list of vertices of a (polyhedral) inverse linear map.

### Input

- `ilm`   -- inverse linear map
- `prune` -- (optional, default: `true`) if `true` removes redundant vertices

### Output

A list of vertices.

### Algorithm

We assume that the underlying set `X` is polyhedral.
Then the result is just the inverse linear map applied to the vertices of `X`.
"""
function vertices_list(ilm::InverseLinearMap; prune::Bool=true)
    # collect vertices list of wrapped set
    vlist_X = vertices_list(ilm.X)

    # create resulting vertices list
    vlist = Vector{eltype(vlist_X)}()
    sizehint!(vlist, length(vlist_X))
    for v in vlist_X
        push!(vlist, ilm.M \ v)
    end

    return prune ? convex_hull(vlist) : vlist
end

"""
    constraints_list(ilm::InverseLinearMap)

Return the list of constraints of a (polyhedral) inverse linear map.

### Input

- `ilm` -- inverse linear map

### Output

The list of constraints of the inverse linear map.

### Notes

We assume that the underlying set `X` is polyhedral, i.e., offers a method
`constraints_list(X)`.

### Algorithm

We fall back to a concrete set representation and apply `linear_map_inverse`.
"""
function constraints_list(ilm::InverseLinearMap)
    return constraints_list(linear_map_inverse(ilm.M, ilm.X))
end

"""
    linear_map(M::AbstractMatrix, ilm::InverseLinearMap)

Return the linear map of a lazy inverse linear map.

### Input

- `M`   -- matrix
- `ilm` -- inverse linear map

### Output

The set representing the linear map of the lazy inverse linear map of a set.

### Notes

This implementation is inefficient because it will compute the
concrete inverse of ``M``, which is what `InverseLinearMap` is supposed to avoid.
"""
function linear_map(M::AbstractMatrix, ilm::InverseLinearMap)
    return linear_map(M * inv(ilm.M), ilm.X)
end

function concretize(ilm::InverseLinearMap)
    return linear_map(inv(ilm.M), concretize(ilm.X))
end
