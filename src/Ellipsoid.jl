export Ellipsoid

"""
    Ellipsoid{N<:Real} <:  AbstractPointSymmetric{N}

Type that represents an ellipsoid.

It is defined as the set

```math
E = \\left\\{ x ∈ \\mathbb{R}^n : (x-c)Q^{-1}(x-c) ≤ 1 \\right\\},
```
where ``c \\in \\mathbb{R}^n`` is its *center* and ``Q \\in \\mathbb{R}^{n×n}``
its *shape matrix*, which should be a positive definite matrix.
An ellipsoid can also be characterized as the image of a Euclidean ball by an
invertible linear transformation.

### Fields

- `center`       -- center of the ellipsoid
- `shape matrix` -- real positive definite matrix, i.e. it is equal to its transpose
                    and ``x^\\mathrm{T}Qx > 0`` for all nonzero ``x``

### Examples

If the center is not specified, it is assumed that the center is the origin.
For instance, a 3D ellipsoid with center at the origin and the shape matrix being
the identity can be created with:

```jldoctest ellipsoid_constructor
julia> E = Ellipsoid(eye(3))
LazySets.Ellipsoid{Float64}([0.0, 0.0, 0.0], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])

julia> dim(E)
3

julia> an_element(E)
3-element Array{Float64,1}:
 0.0
 0.0
 0.0
```

This ellipsoid corresponds to the unit Euclidean ball. Let's evaluate its support
vector in a given direction:

```jldoctest ellipsoid_constructor
julia> σ(ones(3), E)
3-element Array{Float64,1}:
 0.57735
 0.57735
 0.57735
```

A two-dimensional ellipsoid with given center and shape matrix:

```julia
julia> E = Ellipsoid(ones(2), diagm([2.0, 0.5]))
LazySets.Ellipsoid{Float64}([1.0, 1.0], [2.0 0.0; 0.0 0.5])
```
"""
struct Ellipsoid{N<:Real} <: AbstractPointSymmetric{N}
    center::AbstractVector{N}
    shape_matrix::AbstractMatrix{N}

    # default constructor with dimension check
    function Ellipsoid{N}(c::AbstractVector{N}, Q::AbstractMatrix{N}) where {N<:Real}
        @assert length(c) == Base.LinAlg.checksquare(Q)
        return new(c, Q)
    end
end

# type-less convenience constructor
Ellipsoid(c::AbstractVector{N}, Q::AbstractMatrix{N}) where {N<:Real} =
    Ellipsoid{N}(c, Q)

# convenience constructor for ellipsoid centered in the origin
Ellipsoid(Q::AbstractMatrix{N}) where {N<:Real} = Ellipsoid(zeros(N, size(Q, 1)), Q)

"""
    center(E::Ellipsoid{N})::Vector{N} where {N<:AbstractFloat}

Return the center of the ellipsoid.

### Input

- `E` -- ellipsoid

### Output

The center of the ellipsoid.
"""
function center(E::Ellipsoid{N})::Vector{N} where {N<:AbstractFloat}
    return E.center
end

"""
    σ(d::AbstractVector{N},
               E::Ellipsoid{N})::AbstractVector{<:AbstractFloat} where {N<:AbstractFloat}

Return the support vector of an ellipsoid in a given direction.

### Input

- `d` -- direction
- `E` -- ellipsoid

### Output

Support vector in the given direction.

### Algorithm

Let ``E`` be an ellipsoid of center ``c`` and shape matrix ``Q = BB^\\mathrm{T}``.
Its support vector along direction ``d`` can be deduced from that of the unit
Euclidean ball ``\\mathcal{B}_2`` using the algebraic relations for the support
vector,

```math
σ_{B\\mathcal{B}_2 ⊕ c}(d) = c + Bσ_{\\mathcal{B}_2} (B^\\mathrm{T} d)
= c + \\dfrac{Qd}{\\sqrt{d^\\mathrm{T}Q d}}.
```
"""
function σ(d::AbstractVector{N},
           E::Ellipsoid{N})::AbstractVector{<:AbstractFloat} where {N<:AbstractFloat}
    if norm(d, 2) == zero(N)
        return an_element(E)
    end
    α = sqrt(dot(d, E.shape_matrix * d))
    return E.center .+ E.shape_matrix * d ./ α
end

"""
    ∈(x::AbstractVector{N}, E::Ellipsoid{N})::Bool where {N<:AbstractFloat}

Check whether a given point is contained in an ellipsoid.

### Input

- `x` -- point/vector
- `E` -- ellipsoid

### Output

`true` iff `x ∈ E`.
"""
function ∈(x::AbstractVector{N}, E::Ellipsoid{N})::Bool where {N<:AbstractFloat}
    @assert length(x) == dim(E)
    c, Q = E.center, E.shape_matrix
    return dot(x-c, Q \ (x-c)) ≤ 1
end
