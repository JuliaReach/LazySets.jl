import Base: rand,
             ∈,
             isempty

export Line,
       an_element

"""
    Line{N<:Real, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a line in 2D of the form ``a⋅x = b`` (i.e., a special case
of a `Hyperplane`).

### Fields

- `a` -- normal direction (non-zero)
- `b` -- constraint

### Examples

The line ``y = -x + 1``:

```jldoctest
julia> Line([1., 1.], 1.)
Line{Float64,Array{Float64,1}}([1.0, 1.0], 1.0)
```
"""
struct Line{N<:Real, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    a::VN
    b::N

    # default constructor with length constraint
    function Line{N, VN}(a::VN, b::N) where {N<:Real, VN<:AbstractVector{N}}
        @assert length(a) == 2 "lines must be two-dimensional"
        @assert !iszero(a) "a line needs a non-zero normal vector"
        return new{N, VN}(a, b)
    end
end

isoperationtype(::Type{<:Line}) = false

# convenience constructor without type parameter
Line(a::VN, b::N) where {N<:Real, VN<:AbstractVector{N}} = Line{N, VN}(a, b)

# constructor from a LinearConstraint
Line(c::LinearConstraint{N}) where {N<:Real} = Line(c.a, c.b)


"""
    Line(p::AbstractVector{N}, q::AbstractVector{N}) where {N<:Real}

Constructor a line give two points.

### Input

- `p` -- point in 2D
- `q` -- another point in 2D

### Output

The line which passes through `p` and `q`.

### Algorithm

Given two points ``p = (x₁, y₁)`` and ``q = (x₂, y₂)`` the line that passes through these
two points is

```math
ℓ:~~y - y₁ = \\dfrac{(y₂ - y₁)}{(x₂ - x₁)} ⋅ (x-x₁).
```
The particular case ``x₂ = x₁`` defines a line parallel to the ``y``-axis (vertical line).
"""
function Line(p::AbstractVector{N}, q::AbstractVector{N}) where {N<:Real}
    x₁, y₁ = p[1], p[2]
    x₂, y₂ = q[1], q[2]

    if x₁ == x₂  # line is vertical
        a = [one(N), zero(N)]
        b = x₁
        return LazySets.Line(a, b)
    end

    k = (y₁ - y₂)/(x₂ - x₁)
    a = [k, one(N)]
    b = y₁ + k * x₁
    return Line(a, b)
end


# --- polyhedron interface functions ---


"""
    constraints_list(L::Line{N}) where {N<:Real}

Return the list of constraints of a line.

### Input

- `L` -- line

### Output

A list containing two half-spaces.
"""
function constraints_list(L::Line{N}) where {N<:Real}
    return _constraints_list_hyperplane(L.a, L.b)
end


# --- LazySet interface functions ---


"""
    dim(L::Line)::Int

Return the ambient dimension of a line.

### Input

- `L` -- line

### Output

The ambient dimension of the line, which is 2.
"""
function dim(L::Line)::Int
    return 2
end

"""
    σ(d::AbstractVector{N}, L::Line{N}) where {N<:Real}

Return the support vector of a line in a given direction.

### Input

- `d` -- direction
- `L` -- line

### Output

The support vector in the given direction, which is defined the same way as for
the more general `Hyperplane`.
"""
function σ(d::AbstractVector{N}, L::Line{N}) where {N<:Real}
    return σ(d, Hyperplane(L.a, L.b))
end

"""
    isbounded(L::Line)::Bool

Determine whether a line is bounded.

### Input

- `L` -- line

### Output

`false`.
"""
function isbounded(::Line)::Bool
    return false
end

"""
    isuniversal(L::Line{N}, [witness]::Bool=false
               )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}

Check whether a line is universal.

### Input

- `P`       -- line
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ P``

### Algorithm

Witness production falls back to `isuniversal(::Hyperplane)`.
"""
function isuniversal(L::Line{N}, witness::Bool=false
                    )::Union{Bool, Tuple{Bool, Vector{N}}} where {N<:Real}
    if witness
        return isuniversal(Hyperplane(L.a, L.b), true)
    else
        return false
    end
end

"""
    an_element(L::Line{N})::Vector{N} where {N<:Real}

Return some element of a line.

### Input

- `L` -- line

### Output

An element on the line.

### Algorithm

If the ``b`` value of the line is zero, the result is the origin.
Otherwise the result is some ``x = [x1, x2]`` such that ``a·[x1, x2] = b``.
We first find out in which dimension ``a`` is nonzero, say, dimension 1, and
then choose ``x1 = 1`` and accordingly ``x2 = \\frac{b - a1}{a2}``.
"""
function an_element(L::Line{N})::Vector{N} where {N<:Real}
    if L.b == zero(N)
        return zeros(N, 2)
    end
    i = L.a[1] == zero(N) ? 2 : 1
    x = Vector{N}(undef, 2)
    x[3-i] = one(N)
    x[i] = (L.b - L.a[3-i]) / L.a[i]
    return x
end

"""
    ∈(x::AbstractVector{N}, L::Line{N})::Bool where {N<:Real}

Check whether a given point is contained in a line.

### Input

- `x` -- point/vector
- `L` -- line

### Output

`true` iff `x ∈ L`.

### Algorithm

The point ``x`` belongs to the line if and only if ``a⋅x = b`` holds.
"""
function ∈(x::AbstractVector{N}, L::Line{N})::Bool where {N<:Real}
    @assert length(x) == dim(L)
    return dot(L.a, x) == L.b
end

"""
    rand(::Type{Line}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::Line{N}

Create a random line.

### Input

- `Line` -- type for dispatch
- `N`    -- (optional, default: `Float64`) numeric type
- `dim`  -- (optional, default: 2) dimension
- `rng`  -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed` -- (optional, default: `nothing`) seed for reseeding

### Output

A random line.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
Additionally, the constraint `a` is nonzero.
"""
function rand(::Type{Line};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing
             )::Line{N}
    @assert dim == 2 "cannot create a random Line of dimension $dim"
    rng = reseed(rng, seed)
    a = randn(rng, N, dim)
    while iszero(a)
        a = randn(rng, N, dim)
    end
    b = randn(rng, N)
    return Line(a, b)
end

"""
    isempty(L::Line)::Bool

Return if a line is empty or not.

### Input

- `L` -- line

### Output

`false`.
"""
function isempty(L::Line)::Bool
    return false
end

"""
    constrained_dimensions(L::Line{N})::Vector{Int} where {N<:Real}

Return the indices in which a line is constrained.

### Input

- `L` -- line

### Output

A vector of ascending indices `i` such that the line is constrained in dimension
`i`.

### Examples

A line with constraint ``x1 = 0`` is constrained in dimension 1 only.
"""
function constrained_dimensions(L::Line{N})::Vector{Int} where {N<:Real}
    return nonzero_indices(L.a)
end

function _linear_map_hrep(M::AbstractMatrix{N}, P::Line{N}, use_inv::Bool) where {N<:Real}
    constraint = _linear_map_hrep_helper(M, P, use_inv)[1]
    return Line(constraint.a, constraint.b)
end

"""
    translate(L::Line{N}, v::AbstractVector{N}; share::Bool=false
             ) where {N<:Real}

Translate (i.e., shift) a line by a given vector.

### Input

- `L`     -- line
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated line.

### Notes

The normal vector of the line (vector `a` in `a⋅x = b`) is shared with the
original line if `share == true`.

### Algorithm

A line ``a⋅x = b`` is transformed to the line ``a⋅x = b + a⋅v``.
In other words, we add the dot product ``a⋅v`` to ``b``.
"""
function translate(L::Line{N}, v::AbstractVector{N}; share::Bool=false
                  ) where {N<:Real}
    @assert length(v) == dim(L) "cannot translate a $(dim(L))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    a = share ? L.a : copy(L.a)
    b = L.b + dot(L.a, v)
    return Line(a, b)
end
