import Base: rand,
             ∈,
             isempty,
             convert

export HalfSpace, LinearConstraint,
       an_element,
       constrained_dimensions,
       halfspace_left, halfspace_right

"""
    HalfSpace{N<:Real, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a (closed) half-space of the form ``a⋅x ≤ b``.

### Fields

- `a` -- normal direction (non-zero)
- `b` -- constraint

### Examples

The set ``y ≥ 0`` in the plane:

```jldoctest
julia> HalfSpace([0, -1.], 0.)
HalfSpace{Float64,Array{Float64,1}}([0.0, -1.0], 0.0)
```
"""
struct HalfSpace{N<:Real, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    a::VN
    b::N

    function HalfSpace{N, VN}(a::VN, b::N
                             ) where {N<:Real, VN<:AbstractVector{N}}
        @assert !iszero(a) "a half-space needs a non-zero normal vector"
        return new{N, VN}(a, b)
    end
end

isoperationtype(::Type{<:HalfSpace}) = false

# convenience constructor without type parameter
HalfSpace(a::VN, b::N) where {N<:Real, VN<:AbstractVector{N}} =
    HalfSpace{N, VN}(a, b)

"""
    LinearConstraint

Alias for `HalfSpace`
"""
const LinearConstraint = HalfSpace

# --- LazySet interface functions ---


"""
    dim(hs::HalfSpace)::Int

Return the dimension of a half-space.

### Input

- `hs` -- half-space

### Output

The ambient dimension of the half-space.
"""
function dim(hs::HalfSpace)::Int
    return length(hs.a)
end

"""
    ρ(d::AbstractVector{N}, hs::HalfSpace{N})::N where {N<:Real}

Evaluate the support function of a half-space in a given direction.

### Input

- `d`  -- direction
- `hs` -- half-space

### Output

The support function of the half-space.
If the set is unbounded in the given direction, the result is `Inf`.
"""
function ρ(d::AbstractVector{N}, hs::HalfSpace{N})::N where {N<:Real}
    v, unbounded = σ_helper(d, Hyperplane(hs.a, hs.b); error_unbounded=false,
                            halfspace=true)
    if unbounded
        return N(Inf)
    end
    return dot(d, v)
end

"""
    σ(d::AbstractVector{N}, hs::HalfSpace{N}) where {N<:Real}

Return the support vector of a half-space.

### Input

- `d`  -- direction
- `hs` -- half-space

### Output

The support vector in the given direction, which is only defined in the
following two cases:
1. The direction has norm zero.
2. The direction is the half-space's normal direction.
In both cases the result is any point on the boundary (the defining hyperplane).
Otherwise this function throws an error.
"""
function σ(d::AbstractVector{N}, hs::HalfSpace{N}) where {N<:Real}
    v, unbounded = σ_helper(d, Hyperplane(hs.a, hs.b); error_unbounded=true,
                            halfspace=true)
    return v
end

"""
    isbounded(hs::HalfSpace)::Bool

Determine whether a half-space is bounded.

### Input

- `hs` -- half-space

### Output

`false`.
"""
function isbounded(::HalfSpace)::Bool
    return false
end

"""
    isuniversal(hs::HalfSpace{N}, [witness]::Bool=false
               )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a half-space is universal.

### Input

- `P`       -- half-space
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ P``

### Algorithm

Witness production falls back to `isuniversal(::Hyperplane)`.
"""
function isuniversal(hs::HalfSpace{N}, witness::Bool=false
                    )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    if witness
        return isuniversal(Hyperplane(hs.a, hs.b), true)
    else
        return false
    end
end

"""
    an_element(hs::HalfSpace{N})::Vector{N} where {N<:Real}

Return some element of a half-space.

### Input

- `hs` -- half-space

### Output

An element on the defining hyperplane.
"""
function an_element(hs::HalfSpace{N})::Vector{N} where {N<:Real}
    return an_element_helper(Hyperplane(hs.a, hs.b))
end

"""
    ∈(x::AbstractVector{N}, hs::HalfSpace{N})::Bool where {N<:Real}

Check whether a given point is contained in a half-space.

### Input

- `x` -- point/vector
- `hs` -- half-space

### Output

`true` iff ``x ∈ hs``.

### Algorithm

We just check if ``x`` satisfies ``a⋅x ≤ b``.
"""
function ∈(x::AbstractVector{N}, hs::HalfSpace{N})::Bool where {N<:Real}
    return dot(x, hs.a) <= hs.b
end

"""
    rand(::Type{HalfSpace}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::HalfSpace{N}

Create a random half-space.

### Input

- `HalfSpace` -- type for dispatch
- `N`         -- (optional, default: `Float64`) numeric type
- `dim`       -- (optional, default: 2) dimension
- `rng`       -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`      -- (optional, default: `nothing`) seed for reseeding

### Output

A random half-space.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
Additionally, the constraint `a` is nonzero.
"""
function rand(::Type{HalfSpace};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing
             )::HalfSpace{N}
    rng = reseed(rng, seed)
    a = randn(rng, N, dim)
    while iszero(a)
        a = randn(rng, N, dim)
    end
    b = randn(rng, N)
    return HalfSpace(a, b)
end

"""
    isempty(hs::HalfSpace)::Bool

Return if a half-space is empty or not.

### Input

- `hs` -- half-space

### Output

`false`.
"""
function isempty(hs::HalfSpace)::Bool
    return false
end

"""
    constraints_list(hs::HalfSpace{N}) where {N<:Real}

Return the list of constraints of a half-space.

### Input

- `hs` -- half-space

### Output

A singleton list containing the half-space.
"""
function constraints_list(hs::HalfSpace{N}) where {N<:Real}
    return [hs]
end

"""
    constraints_list(A::AbstractMatrix{N}, b::AbstractVector{N}) where {N<:Real}

Convert a matrix-vector representation to a linear-constraint representation.

### Input

- `A` -- matrix
- `b` -- vector

### Output

A list of linear constraints.
"""
function constraints_list(A::AbstractMatrix{N}, b::AbstractVector{N}
                         ) where {N<:Real}
    m = size(A, 1)
    @assert m == length(b) "a matrix with $m rows is incompatible with a " *
                           "vector of length $(length(b))"

    return [LinearConstraint(A[i, :], b[i]) for i in 1:m]
end

"""
    constrained_dimensions(hs::HalfSpace{N})::Vector{Int} where {N<:Real}

Return the indices in which a half-space is constrained.

### Input

- `hs` -- half-space

### Output

A vector of ascending indices `i` such that the half-space is constrained in
dimension `i`.

### Examples

A 2D half-space with constraint ``x1 ≥ 0`` is constrained in dimension 1 only.
"""
function constrained_dimensions(hs::HalfSpace{N})::Vector{Int} where {N<:Real}
    return nonzero_indices(hs.a)
end

"""
    halfspace_left(p::AbstractVector{N},
                   q::AbstractVector{N})::HalfSpace{N} where {N<:Real}

Return a half-space describing the 'left' of a two-dimensional line segment
through two points.

### Input

- `p` -- first point
- `q` -- second point

### Output

The half-space whose boundary goes through the two points `p` and `q` and which
describes the left-hand side of the directed line segment `pq`.

### Algorithm

The implementation is simple: the half-space ``a⋅x ≤ b`` is calculated as
`a = [dy, -dx]`, where ``d = (dx, dy)`` denotes the line segment
`pq`, that is, ``\\vec{d} = \\vec{p} - \\vec{q}``, and `b = dot(p, a)`.

### Examples

The left half-space of the "east" and "west" directions in two-dimensions are
the upper and lower half-spaces:

```jldoctest halfspace_left
julia> using LazySets: halfspace_left

julia> halfspace_left([0.0, 0.0], [1.0, 0.0])
HalfSpace{Float64,Array{Float64,1}}([0.0, -1.0], 0.0)

julia> halfspace_left([0.0, 0.0], [-1.0, 0.0])
HalfSpace{Float64,Array{Float64,1}}([0.0, 1.0], 0.0)
```

We create a box from the sequence of line segments that describe its edges:

```jldoctest halfspace_left
julia> H1 = halfspace_left([-1.0, -1.0], [1.0, -1.0]);

julia> H2 = halfspace_left([1.0, -1.0], [1.0, 1.0]);

julia> H3 = halfspace_left([1.0, 1.0], [-1.0, 1.0]);

julia> H4 = halfspace_left([-1.0, 1.0], [-1.0, -1.0]);

julia> H = HPolygon([H1, H2, H3, H4]);

julia> B = BallInf([0.0, 0.0], 1.0);

julia> B ⊆ H && H ⊆ B
true
```
"""
function halfspace_left(p::AbstractVector{N},
                        q::AbstractVector{N})::HalfSpace{N} where {N<:Real}
    @assert length(p) == length(q) == 2 "the points must be two-dimensional"
    @assert p != q "the points must not be equal"
    a = [q[2] - p[2], p[1] - q[1]]
    return HalfSpace(a, dot(p, a))
end

"""
    halfspace_right(p::AbstractVector{N},
                    q::AbstractVector{N})::HalfSpace{N} where {N<:Real}

Return a half-space describing the 'right' of a two-dimensional line segment
through two points.

### Input

- `p` -- first point
- `q` -- second point

### Output

The half-space whose boundary goes through the two points `p` and `q` and which
describes the right-hand side of the directed line segment `pq`.

### Algorithm

See the documentation of `halfspace_left`.
"""
function halfspace_right(p::AbstractVector{N},
                         q::AbstractVector{N})::HalfSpace{N} where {N<:Real}
    return halfspace_left(q, p)
end

"""
    is_tighter_same_dir_2D(c1::LinearConstraint{N},
                           c2::LinearConstraint{N}) where {N<:Real}

Check if the first of two two-dimensional constraints with equivalent normal
direction is tighter.

### Input

- `c1`     -- first linear constraint
- `c2`     -- second linear constraint
- `strict` -- (optional; default: `false`) check for strictly tighter
              constraints?

### Output

`true` iff the first constraint is tighter.
"""
function is_tighter_same_dir_2D(c1::LinearConstraint{N},
                                c2::LinearConstraint{N};
                                strict::Bool=false) where {N<:Real}
    @assert dim(c1) == dim(c2) == 2 "this method requires 2D constraints"
    @assert c1.a <= c2.a <= c1.a "the constraints must have the same " *
        "normal direction"

    lt = strict ? (<) : (<=)
    if c1.a[1] == zero(N)
        @assert c2.a[1] == zero(N)
        return lt(c1.b, c1.a[2] / c2.a[2] * c2.b)
    end
    return lt(c1.b, c1.a[1] / c2.a[1] * c2.b)
end

function _linear_map_hrep(M::AbstractMatrix{N}, P::HalfSpace{N}, use_inv::Bool) where {N<:Real}
    constraint = _linear_map_hrep_helper(M, P, use_inv)[1]
    return HalfSpace(constraint.a, constraint.b)
end

"""
    translate(hs::HalfSpace{N}, v::AbstractVector{N}; share::Bool=false
             ) where {N<:Real}

Translate (i.e., shift) a half-space by a given vector.

### Input

- `hs`    -- half-space
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated half-space.

### Notes

The normal vectors of the halfspace (vector `a` in `a⋅x ≤ b`) is shared with the
original halfspace if `share == true`.

### Algorithm

A half-space ``a⋅x ≤ b`` is transformed to the half-space ``a⋅x ≤ b + a⋅v``.
In other words, we add the dot product ``a⋅v`` to ``b``.
"""
function translate(hs::HalfSpace{N}, v::AbstractVector{N}; share::Bool=false
                  ) where {N<:Real}
    @assert length(v) == dim(hs) "cannot translate a $(dim(hs))-dimensional " *
                                 "set by a $(length(v))-dimensional vector"
    a = share ? hs.a : copy(hs.a)
    b = hs.b + dot(hs.a, v)
    return HalfSpace(a, b)
end
