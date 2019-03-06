import Base.rand

export Ellipsoid

"""
    Ellipsoid{N<:AbstractFloat} <:  AbstractCentrallySymmetric{N}

Type that represents an ellipsoid.

It is defined as the set

```math
E = \\left\\{ x ∈ \\mathbb{R}^n : (x-c)Q^{-1}(x-c) ≤ 1 \\right\\},
```
where ``c \\in \\mathbb{R}^n`` is its *center* and ``Q \\in \\mathbb{R}^{n×n}``
its *shape matrix*, which should be a positive definite matrix.
An ellipsoid can also be characterized as the image of a Euclidean ball by an
invertible linear transformation. It is the higher-dimensional generalization
of an ellipse.

### Fields

- `center`       -- center of the ellipsoid
- `shape matrix` -- real positive definite matrix, i.e. it is equal to its transpose
                    and ``x^\\mathrm{T}Qx > 0`` for all nonzero ``x``

### Examples

If the center is not specified, it is assumed that the center is the origin.
For instance, a 3D ellipsoid with center at the origin and the shape matrix being
the identity can be created with:

```jldoctest ellipsoid_constructor
julia> E = Ellipsoid(Matrix{Float64}(I, 3, 3))
Ellipsoid{Float64}([0.0, 0.0, 0.0], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])

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
 0.5773502691896258
 0.5773502691896258
 0.5773502691896258
```

A two-dimensional ellipsoid with given center and shape matrix:

```julia
julia> E = Ellipsoid(ones(2), Diagonal([2.0, 0.5]))
Ellipsoid{Float64}([1.0, 1.0], [2.0 0.0; 0.0 0.5])
```
"""
struct Ellipsoid{N<:AbstractFloat} <: AbstractCentrallySymmetric{N}
    center::AbstractVector{N}
    shape_matrix::AbstractMatrix{N}

    # default constructor with dimension check
    function Ellipsoid{N}(c::AbstractVector{N},
                          Q::AbstractMatrix{N}) where {N<:AbstractFloat}
        @assert length(c) == checksquare(Q)
        @assert isposdef(Q) "an ellipsoid's shape matrix must be positive " *
            "definite"
        return new{N}(c, Q)
    end
end

# convenience constructor without type parameter
Ellipsoid(c::AbstractVector{N}, Q::AbstractMatrix{N}) where {N<:AbstractFloat} =
    Ellipsoid{N}(c, Q)

# convenience constructor for ellipsoid centered in the origin
Ellipsoid(Q::AbstractMatrix{N}) where {N<:AbstractFloat} =
    Ellipsoid(zeros(N, size(Q, 1)), Q)


# --- AbstractCentrallySymmetric interface functions ---


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


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, E::Ellipsoid{N}) where {N<:AbstractFloat}

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
function σ(d::AbstractVector{N}, E::Ellipsoid{N}) where {N<:AbstractFloat}
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

### Algorithm

The point ``x`` belongs to the ellipsoid of center ``c`` and shape matrix ``Q``
if and only if

```math
(x-c)^\\mathrm{T} Q^{-1} (x-c) ≤ 1.
```
"""
function ∈(x::AbstractVector{N}, E::Ellipsoid{N})::Bool where {N<:AbstractFloat}
    @assert length(x) == dim(E)
    w, Q = x-E.center, E.shape_matrix
    return dot(w, Q \ w) ≤ 1
end

"""
    rand(::Type{Ellipsoid}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::Ellipsoid{N}

Create a random ellipsoid.

### Input

- `Ellipsoid` -- type for dispatch
- `N`         -- (optional, default: `Float64`) numeric type
- `dim`       -- (optional, default: 2) dimension
- `rng`       -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`      -- (optional, default: `nothing`) seed for reseeding

### Output

A random ellipsoid.

### Algorithm

The center is a normally distributed vector with entries of mean 0 and standard
deviation 1.

The idea for the shape matrix comes from
[here](https://math.stackexchange.com/a/358092).
The matrix is symmetric positive definite, but also diagonally dominant.

```math
Q =  \frac{1}{2}(S + S^T) + nI,
```
where ``n`` = `dim` (defaults to 2), and ``S`` is a ``n \\times n`` random
matrix whose coefficients are uniformly distributed in the interval ``[-1, 1]``.
"""
function rand(::Type{Ellipsoid};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing
             )::Ellipsoid{N}
    rng = reseed(rng, seed)
    center = randn(rng, N, dim)
    # random entries in [-1, 1]
    # this needs a bit of code because 'rand' only samples from [0, 1]
    shape_matrix = Matrix{N}(undef, dim, dim)
    for j in 1:dim
        for i in 1:dim
            entry = rand(N)
            if rand(Bool)
                entry = -entry
            end
            shape_matrix[i, j] = entry
        end
    end
    # make diagonally dominant
    shape_matrix = N(0.5) * (shape_matrix + shape_matrix') +
                   Matrix{N}(dim*I, dim, dim)
    return Ellipsoid(center, shape_matrix)
end

"""
    translate(E::Ellipsoid{N}, v::AbstractVector{N}; share::Bool=false
             ) where {N<:AbstractFloat}

Translate (i.e., shift) an ellipsoid by a given vector.

### Input

- `E`     -- ellipsoid
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated ellipsoid.

### Notes

The shape matrix is shared with the original ellipsoid if `share == true`.

### Algorithm

We add the vector to the center of the ellipsoid.
"""
function translate(E::Ellipsoid{N}, v::AbstractVector{N}; share::Bool=false
                  ) where {N<:AbstractFloat}
    @assert length(v) == dim(E) "cannot translate a $(dim(E))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = center(E) + v
    shape_matrix = share ? E.shape_matrix : copy(E.shape_matrix)
    return Ellipsoid(c, shape_matrix)
end
