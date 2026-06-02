import Base: *

"""
    InverseLinearMap{N, S<:LazySet{N}, NM, MAT<:AbstractMatrix{NM}}
        <: AbstractAffineMap{N, S}

Given a linear transformation ``M``, this type represents the linear
transformation ``M⁻¹⋅X`` of a set ``X`` without actually computing ``M⁻¹``.

### Fields

- `M` -- matrix (typically invertible, which can be checked in the constructor)
- `X` -- set

### Notes

Many set operations avoid computing the inverse of the matrix.

In principle, the matrix does not have to be invertible (it can for instance be
rectangular) for many set operations.

This type is parametric in the elements of the inverse linear map, `NM`, which
is independent of the numeric type of the wrapped set (`N`). Typically `NM = N`,
but there may be exceptions, e.g., if `NM` is an interval that holds numbers
of type `N`, where `N` is a floating-point type such as `Float64`.

### Examples

For the examples we create a ``3×3`` matrix and a unit three-dimensional square.

```jldoctest ilp_constructor
julia> A = [1 2 3; 2 3 1; 3 1 2];

julia> X = BallInf([0, 0, 0], 1);

julia> ilm = InverseLinearMap(A, X)
InverseLinearMap{Int64, BallInf{Int64, Vector{Int64}}, Int64, Matrix{Int64}}([1 2 3; 2 3 1; 3 1 2], BallInf{Int64, Vector{Int64}}([0, 0, 0], 1))
```

Applying an inverse linear map to a `InverseLinearMap` object combines the two maps into
a single `InverseLinearMap` instance.

```jldoctest ilp_constructor
julia> B = transpose(A); ilm2 = InverseLinearMap(B, ilm)
InverseLinearMap{Int64, BallInf{Int64, Vector{Int64}}, Int64, Matrix{Int64}}([14 11 11; 11 14 11; 11 11 14], BallInf{Int64, Vector{Int64}}([0, 0, 0], 1))

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
struct InverseLinearMap{N,S<:LazySet{N},NM,MAT<:AbstractMatrix{NM}} <: AbstractAffineMap{N,S}
    M::MAT
    X::S

    # default constructor with dimension match check
    function InverseLinearMap(M::MAT, X::S;
                              check_invertibility::Bool=false) where {N,S<:LazySet{N},NM,
                                                                      MAT<:AbstractMatrix{NM}}
        @assert dim(X) == size(M, 1) "a linear map of size $(size(M)) cannot " *
                                     "be applied to a set of dimension $(dim(X))"
        if check_invertibility
            @assert isinvertible(M) "the linear map is not invertible"
        end
        return new{N,S,NM,MAT}(M, X)
    end
end

# convenience constructor from a UniformScaling
function InverseLinearMap(M::UniformScaling, X::LazySet;
                          check_invertibility::Bool=false)
    return InverseLinearMap(M.λ, X; check_invertibility=check_invertibility)
end

# convenience constructor from a scalar
function InverseLinearMap(α::Real, X::LazySet; check_invertibility::Bool=false)
    if check_invertibility
        @assert !iszero(α) "the linear map is not invertible"
    end

    if isone(α)
        return X
    end

    D = Diagonal(fill(α, dim(X)))
    return InverseLinearMap(D, X; check_invertibility=false)
end

# combine two inverse linear maps into a single inverse linear map
function InverseLinearMap(M::AbstractMatrix, ilm::InverseLinearMap)
    return InverseLinearMap(ilm.M * M, ilm.X)
end

# ZeroSet is almost absorbing for InverseLinearMap (only the dimension changes)
function InverseLinearMap(M::AbstractMatrix, Z::ZeroSet{N}) where {N}
    @assert dim(Z) == size(M, 2) "a linear map of size $(size(M)) cannot " *
                                 "be applied to a set of dimension $(dim(Z))"
    return ZeroSet{N}(size(M, 1))
end

# EmptySet is almost absorbing for InverseLinearMap (only the dimension changes)
function InverseLinearMap(M::AbstractMatrix, ∅::EmptySet{N}) where {N}
    @assert dim(∅) == size(M, 2) "a linear map of size $(size(M)) cannot " *
                                 "be applied to a set of dimension $(dim(∅))"
    return EmptySet{N}(size(M, 1))
end

function matrix(ilm::InverseLinearMap)
    return inv(ilm.M)
end

function vector(ilm::InverseLinearMap{N}) where {N}
    return spzeros(N, dim(ilm))
end

function set(ilm::InverseLinearMap)
    return ilm.X
end

include("an_element.jl")
include("concretize.jl")
include("constraints_list.jl")
include("dim.jl")
include("vertices_list.jl")
include("in.jl")
include("linear_map.jl")
include("support_function.jl")
include("support_vector.jl")
