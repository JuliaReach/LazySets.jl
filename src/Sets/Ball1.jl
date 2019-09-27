import Base: rand,
             ∈

export Ball1

"""
    Ball1{N<:Real} <: AbstractCentrallySymmetricPolytope{N}

Type that represents a ball in the 1-norm (also known as the Manhattan norm).
The ball is also known as a
[cross-polytope](https://en.wikipedia.org/wiki/Cross-polytope).

It is defined as the set

```math
\\mathcal{B}_1^n(c, r) = \\{ x ∈ \\mathbb{R}^n : ∑_{i=1}^n |c_i - x_i| ≤ r \\},
```
where ``c ∈ \\mathbb{R}^n`` is its center and ``r ∈ \\mathbb{R}_+`` its radius.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a scalar (``≥ 0``)

### Examples

Unit ball in the 1-norm in the plane:

```jldoctest ball1_constructor
julia> B = Ball1(zeros(2), 1.)
Ball1{Float64}([0.0, 0.0], 1.0)
julia> dim(B)
2
```

We evaluate the support vector in the East direction:

```jldoctest ball1_constructor
julia> σ([0.,1], B)
2-element Array{Float64,1}:
 0.0
 1.0
```
"""
struct Ball1{N<:Real} <: AbstractCentrallySymmetricPolytope{N}
    center::Vector{N}
    radius::N

    # default constructor with domain constraint for radius
    function Ball1{N}(center::Vector{N}, radius::N) where {N<:Real}
        @assert radius >= zero(N) "radius must not be negative"
        return new{N}(center, radius)
    end
end

isoperationtype(::Type{<:Ball1}) = false

# convenience constructor without type parameter
Ball1(center::Vector{N}, radius::N) where {N<:Real} = Ball1{N}(center, radius)


# --- AbstractCentrallySymmetric interface functions ---


"""
    center(B::Ball1{N})::Vector{N} where {N<:Real}

Return the center of a ball in the 1-norm.

### Input

- `B` -- ball in the 1-norm

### Output

The center of the ball in the 1-norm.
"""
function center(B::Ball1{N})::Vector{N} where {N<:Real}
    return B.center
end


# --- AbstractPolytope interface functions ---


"""
    vertices_list(B::Ball1{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a ball in the 1-norm.

### Input

- `B` -- ball in the 1-norm

### Output

A list containing the vertices of the ball in the 1-norm.
"""
function vertices_list(B::Ball1{N})::Vector{Vector{N}} where {N<:Real}
    # fast evaluation if B has radius 0
    if iszero(B.radius)
        return [B.center]
    end
    vertices = Vector{Vector{N}}()
    sizehint!(vertices, 2 * dim(B))
    v = copy(B.center)
    for i in 1:dim(B)
        v[i] += B.radius
        push!(vertices, copy(v))
        v[i] = B.center[i] - B.radius
        push!(vertices, copy(v))
        v[i] = B.center[i]
    end
    return vertices
end


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, B::Ball1{N}) where {N<:Real}

Return the support vector of a ball in the 1-norm in a given direction.

### Input

- `d` -- direction
- `B` -- ball in the 1-norm

### Output

Support vector in the given direction.
"""
function σ(d::AbstractVector{N}, B::Ball1{N}) where {N<:Real}
    res = copy(B.center)
    imax = argmax(abs.(d))
    res[imax] += sign(d[imax]) * B.radius
    return res
end

"""
    ∈(x::AbstractVector{N}, B::Ball1{N})::Bool where {N<:Real}

Check whether a given point is contained in a ball in the 1-norm.

### Input

- `x` -- point/vector
- `B` -- ball in the 1-norm

### Output

`true` iff ``x ∈ B``.

### Notes

This implementation is worst-case optimized, i.e., it is optimistic and first
computes (see below) the whole sum before comparing to the radius.
In applications where the point is typically far away from the ball, a fail-fast
implementation with interleaved comparisons could be more efficient.

### Algorithm

Let ``B`` be an ``n``-dimensional ball in the 1-norm with radius ``r`` and let
``c_i`` and ``x_i`` be the ball's center and the vector ``x`` in dimension
``i``, respectively.
Then ``x ∈ B`` iff ``∑_{i=1}^n |c_i - x_i| ≤ r``.

### Examples

```jldoctest
julia> B = Ball1([1., 1.], 1.);

julia> [.5, -.5] ∈ B
false
julia> [.5, 1.5] ∈ B
true
```
"""
function ∈(x::AbstractVector{N}, B::Ball1{N})::Bool where {N<:Real}
    @assert length(x) == dim(B)
    sum = zero(N)
    for i in eachindex(x)
        sum += abs(B.center[i] - x[i])
    end
    return sum <= B.radius
end

"""
    rand(::Type{Ball1}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::Ball1{N}

Create a random ball in the 1-norm.

### Input

- `Ball1` -- type for dispatch
- `N`     -- (optional, default: `Float64`) numeric type
- `dim`   -- (optional, default: 2) dimension
- `rng`   -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`  -- (optional, default: `nothing`) seed for reseeding

### Output

A random ball in the 1-norm.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
Additionally, the radius is nonnegative.
"""
function rand(::Type{Ball1};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing
             )::Ball1{N}
    rng = reseed(rng, seed)
    center = randn(rng, N, dim)
    radius = abs(randn(rng, N))
    return Ball1(center, radius)
end

"""
    constraints_list(P::Ball1{N}) where {N<:Real}

Return the list of constraints defining a ball in the 1-norm.

### Input

- `B` -- ball in the 1-norm

### Output

The list of constraints of the ball.

### Algorithm

The constraints can be defined as ``d_i^T (x-c) ≤ r`` for all ``d_i``, where
``d_i`` is a vector with elements ``1`` or ``-1`` in ``n`` dimensions. To span
all possible ``d_i``, the function `Iterators.product` is used.
"""
function constraints_list(B::Ball1{N}) where {N<:Real}
    n = LazySets.dim(B)
    c, r = B.center, B.radius
    clist = Vector{LinearConstraint{N}}(undef, 2^n)
    for (i, di) in enumerate(Iterators.product([[one(N), -one(N)] for i = 1:n]...))
        d = collect(di) # tuple -> vector
        clist[i] = LinearConstraint(d, dot(d, c) + r)
    end
    return clist
end

"""
    translate(B::Ball1{N}, v::AbstractVector{N}) where {N<:Real}

Translate (i.e., shift) a ball in the 1-norm by a given vector.

### Input

- `B` -- ball in the 1-norm
- `v` -- translation vector

### Output

A translated ball in the 1-norm.

### Algorithm

We add the vector to the center of the ball.
"""
function translate(B::Ball1{N}, v::AbstractVector{N}) where {N<:Real}
    @assert length(v) == dim(B) "cannot translate a $(dim(B))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return Ball1(center(B) + v, B.radius)
end
