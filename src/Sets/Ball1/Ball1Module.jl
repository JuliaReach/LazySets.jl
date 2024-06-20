export Ball1

"""
    Ball1{N, VN<:AbstractVector{N}} <: AbstractCentrallySymmetricPolytope{N}

Type that represents a ball in the 1-norm (also known as the Manhattan norm).
The ball is also known as a
[cross-polytope](https://en.wikipedia.org/wiki/Cross-polytope).

It is defined as the set

```math
\\mathcal{B}_1^n(c, r) = \\{ x ∈ ℝ^n : ∑_{i=1}^n |c_i - x_i| ≤ r \\},
```
where ``c ∈ ℝ^n`` is its center and ``r ∈ ℝ_+`` its radius.

### Fields

- `center` -- center of the ball as a real vector
- `radius` -- radius of the ball as a scalar (``≥ 0``)

### Examples

The unit ball in the 1-norm in the plane:

```jldoctest ball1_constructor
julia> B = Ball1(zeros(2), 1.0)
Ball1{Float64, Vector{Float64}}([0.0, 0.0], 1.0)
julia> dim(B)
2
```

We evaluate the support vector in the North direction:

```jldoctest ball1_constructor
julia> σ([0.0, 1.0], B)
2-element Vector{Float64}:
 0.0
 1.0
```
"""
struct Ball1{N,VN<:AbstractVector{N}} <: AbstractCentrallySymmetricPolytope{N}
    center::VN
    radius::N

    # default constructor with domain constraint for radius
    function Ball1(center::VN, radius::N) where {N,VN<:AbstractVector{N}}
        @assert radius >= zero(N) "the radius must not be negative"
        return new{N,VN}(center, radius)
    end
end

isoperationtype(::Type{<:Ball1}) = false

_vector_type(B::Ball1{N,VN}) where {N,VN} = VN

"""
    center(B::Ball1)

Return the center of a ball in the 1-norm.

### Input

- `B` -- ball in the 1-norm

### Output

The center of the ball in the 1-norm.
"""
function center(B::Ball1)
    return B.center
end

function radius_ball(B::Ball1)
    return B.radius
end

function ball_norm(B::Ball1)
    N = eltype(B)
    return one(N)
end

function low(B::Ball1)
    return _low_AbstractBallp(B)
end

function low(B::Ball1, i::Int)
    return _low_AbstractBallp(B, i)
end

function high(B::Ball1)
    return _high_AbstractBallp(B)
end

function high(B::Ball1, i::Int)
    return _high_AbstractBallp(B, i)
end

"""
    vertices_list(B::Ball1)

Return the list of vertices of a ball in the 1-norm.

### Input

- `B` -- ball in the 1-norm

### Output

A list containing the vertices of the ball in the 1-norm.

### Notes

In ``n`` dimensions there are ``2n`` vertices (unless the radius is 0).
"""
function vertices_list(B::Ball1)
    # fast evaluation if B has radius 0
    if iszero(B.radius)
        return [B.center]
    end
    VN = _vector_type(B)
    vertices = Vector{VN}(undef, 2 * dim(B))
    j = 0
    v = copy(B.center)
    @inbounds for i in 1:dim(B)
        ci = v[i]
        v[i] += B.radius
        j += 1
        vertices[j] = copy(v)
        v[i] = ci - B.radius
        j += 1
        vertices[j] = copy(v)
        v[i] = ci  # restore old value
    end
    return vertices
end

"""
    σ(d::AbstractVector, B::Ball1)

Return the support vector of a ball in the 1-norm in the given direction.

### Input

- `d` -- direction
- `B` -- ball in the 1-norm

### Output

The support vector in the given direction.
"""
function σ(d::AbstractVector, B::Ball1)
    res = copy(B.center)
    imax = argmaxabs(d)
    res[imax] += sign(d[imax]) * B.radius
    return res
end

"""
    ρ(d::AbstractVector, B::Ball1)

Evaluate the support function of a ball in the 1-norm in the given direction.

### Input

- `d` -- direction
- `B` -- ball in the 1-norm

### Output

Evaluation of the support function in the given direction.

### Algorithm

Let ``c`` and ``r`` be the center and radius of the ball ``B`` in the 1-norm,
respectively. Then:

```math
ρ(d, B) = ⟨d, c⟩ + r ‖d‖_∞.
```
"""
function ρ(d::AbstractVector, B::Ball1)
    return dot(d, B.center) + B.radius * maximum(abs, d)
end

"""
    ∈(x::AbstractVector, B::Ball1, [failfast]::Bool=false)

Check whether a given point is contained in a ball in the 1-norm.

### Input

- `x` -- point/vector
- `B` -- ball in the 1-norm
- `failfast` -- (optional, default: `false`) optimization for negative answer

### Output

`true` iff ``x ∈ B``.

### Notes

The default behavior (`failfast == false`) is worst-case optimized, i.e., the
implementation is optimistic and first computes (see below) the whole sum before
comparing to the radius.
In applications where the point is typically far away from the ball, the option
`failfast == true` terminates faster.

### Algorithm

Let ``B`` be an ``n``-dimensional ball in the 1-norm with radius ``r`` and let
``c_i`` and ``x_i`` be the ball's center and the vector ``x`` in dimension
``i``, respectively.
Then ``x ∈ B`` iff ``∑_{i=1}^n |c_i - x_i| ≤ r``.

### Examples

```jldoctest
julia> B = Ball1([1.0, 1.0], 1.0);

julia> [0.5, -0.5] ∈ B
false
julia> [0.5, 1.5] ∈ B
true
```
"""
function ∈(x::AbstractVector, B::Ball1, failfast::Bool=false)
    @assert length(x) == dim(B) "a $(length(x))-dimensional vector is " *
                                "incompatible with a $(dim(B))-dimensional set"
    N = promote_type(eltype(x), eltype(B))
    sum = zero(N)
    @inbounds for (i, xi) in enumerate(x)
        sum += abs(B.center[i] - xi)
        if failfast && sum > B.radius
            return false
        end
    end
    return sum <= B.radius
end

"""
    rand(::Type{Ball1}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )

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
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    center = randn(rng, N, dim)
    radius = abs(randn(rng, N))
    return Ball1(center, radius)
end

"""
    constraints_list(P::Ball1)

Return the list of constraints of a ball in the 1-norm.

### Input

- `B` -- ball in the 1-norm

### Output

The list of constraints of the ball.

### Notes

In ``n`` dimensions there are ``2^n`` constraints (unless the radius is 0).

### Algorithm

The constraints can be defined as ``d_i^T (x-c) ≤ r`` for all ``d_i``, where
``d_i`` is a vector with elements ``1`` or ``-1`` in ``n`` dimensions. To span
all possible ``d_i``, the function `Iterators.product` is used.
"""
function constraints_list(B::Ball1)
    n = dim(B)
    c, r = B.center, B.radius
    N = eltype(B)
    clist = Vector{HalfSpace{N,Vector{N}}}(undef, 2^n)
    for (i, di) in enumerate(Iterators.product([[one(N), -one(N)] for i in 1:n]...))
        d = collect(di) # tuple -> vector
        clist[i] = HalfSpace(d, dot(d, c) + r)
    end
    return clist
end

"""
    translate!(B::Ball1, v::AbstractVector)

Translate (i.e., shift) a ball in the 1-norm by the given vector, in-place.

### Input

- `B` -- ball in the 1-norm
- `v` -- translation vector

### Output

The in-place translated ball in the 1-norm.

### Algorithm

We add the vector to the center of the ball.

### Notes

See also [`translate(::Ball1, ::AbstractVector)`](@ref) for the out-of-place version.
"""
function translate!(B::Ball1, v::AbstractVector)
    @assert length(v) == dim(B) "cannot translate a $(dim(B))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = B.center
    c .+= v
    return B
end

function project(B::Ball1, block::AbstractVector{Int}; kwargs...)
    return Ball1(B.center[block], B.radius)
end

"""
    reflect(B::Ball1)

Concrete reflection of a ball in the 1-norm `B`, resulting in the reflected set
`-B`.

### Input

- `B` -- ball in the 1-norm

### Output

The `Ball1` representing `-B`.

### Algorithm

If ``B`` has center ``c`` and radius ``r``, then ``-B`` has center ``-c`` and
radius ``r``.
"""
function reflect(B::Ball1)
    return Ball1(-center(B), B.radius)
end

function scale(α::Real, B::Ball1)
    return Ball1(B.center .* α, B.radius * abs(α))
end
