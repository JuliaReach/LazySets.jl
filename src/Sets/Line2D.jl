import Base: rand,
             ∈,
             isempty

export Line2D,
       an_element

"""
    Line2D{N, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a line in 2D of the form ``a⋅x = b`` (i.e., a special case
of a `Hyperplane`).

### Fields

- `a` -- normal direction (non-zero)
- `b` -- constraint

### Examples

The line ``y = -x + 1``:

```jldoctest
julia> Line2D([1., 1.], 1.)
Line2D{Float64, Vector{Float64}}([1.0, 1.0], 1.0)
```

The alternative constructor takes two 2D points (`AbstractVector`s) `p` and `q`
and creates a canonical line from `p` to `q`. See the algorithm section below
for details.

```jldoctest
julia> Line2D([1., 1.], [2., 2])
Line2D{Float64, Vector{Float64}}([-1.0, 1.0], 0.0)
```

### Algorithm

Given two points ``p = (x₁, y₁)`` and ``q = (x₂, y₂)``, the line that passes
through these points is

```math
ℓ:~~y - y₁ = \\dfrac{(y₂ - y₁)}{(x₂ - x₁)} ⋅ (x-x₁).
```
The particular case ``x₂ = x₁`` defines a line parallel to the ``y``-axis
(vertical line).
"""
struct Line2D{N,VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    a::VN
    b::N

    # default constructor with length constraint
    function Line2D(a::VN, b::N) where {N,VN<:AbstractVector{N}}
        @assert length(a) == 2 "a Line2D must be two-dimensional"
        @assert !iszero(a) "a line needs a non-zero normal vector"
        return new{N,VN}(a, b)
    end
end

function Line2D(p::AbstractVector, q::AbstractVector)
    @assert length(p) == length(q) == 2 "a Line2D must be two-dimensional"

    N = promote_type(eltype(p), eltype(q))
    x₁, y₁ = @inbounds p[1], p[2]
    x₂, y₂ = @inbounds q[1], q[2]

    if x₁ == x₂  # line is vertical
        @assert y₁ != y₂ "a line needs two distinct points"
        a = [one(N), zero(N)]
        b = x₁
        return Line2D(a, b)
    end

    k = (y₁ - y₂) / (x₂ - x₁)
    a = [k, one(N)]
    b = y₁ + k * x₁
    return Line2D(a, b)
end

isoperationtype(::Type{<:Line2D}) = false

"""
    constraints_list(L::Line2D)

Return the list of constraints of a 2D line.

### Input

- `L` -- 2D line

### Output

A list containing two half-spaces.
"""
function constraints_list(L::Line2D)
    return _constraints_list_hyperplane(L.a, L.b)
end

"""
    dim(L::Line2D)

Return the ambient dimension of a 2D line.

### Input

- `L` -- 2D line

### Output

The ambient dimension of the line, which is 2.
"""
function dim(L::Line2D)
    return 2
end

"""
    σ(d::AbstractVector, L::Line2D)

Return the support vector of a 2D line in a given direction.

### Input

- `d` -- direction
- `L` -- 2D line

### Output

The support vector in the given direction, which is defined the same way as for
the more general `Hyperplane`.
"""
function σ(d::AbstractVector, L::Line2D)
    v, unbounded = _σ_hyperplane_halfspace(d, L.a, L.b; error_unbounded=true,
                                           halfspace=false)
    return v
end

"""
    isbounded(L::Line2D)

Check whether a 2D line is bounded.

### Input

- `L` -- 2D line

### Output

`false`.
"""
function isbounded(::Line2D)
    return false
end

"""
    isuniversal(L::Line2D, [witness]::Bool=false)

Check whether a 2D line is universal.

### Input

- `L`       -- 2D line
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ L``

### Algorithm

Witness production falls back to `isuniversal(::Hyperplane)`.
"""
function isuniversal(L::Line2D, witness::Bool=false)
    if witness
        v = _non_element_halfspace(L.a, L.b)
        return (false, v)
    else
        return false
    end
end

"""
    an_element(L::Line2D)

Return some element of a 2D line.

### Input

- `L` -- 2D line

### Output

An element on the line.

### Algorithm

The algorithm is a 2D specialization of the `Hyperplane` algorithm.

If the ``b`` value of the line is zero, the result is the origin.
Otherwise the result is some ``x = (x_1, x_2)ᵀ`` such that ``a·x = b``.
We first find out the dimension ``i`` in which ``a = (a_1, a_2)ᵀ`` is nonzero
and then choose ``x_i = \\frac{b}{a_i}`` and ``x_{3-i} = 0``.
"""
function an_element(L::Line2D)
    N = eltype(L)
    if iszero(L.b)
        return zeros(N, 2)
    end
    i = @inbounds iszero(L.a[1]) ? 2 : 1
    x = Vector{N}(undef, 2)
    @inbounds x[i] = L.b / L.a[i]
    @inbounds x[3 - i] = zero(N)
    return x
end

"""
    ∈(x::AbstractVector, L::Line2D)

Check whether a given point is contained in a 2D line.

### Input

- `x` -- point/vector
- `L` -- 2D line

### Output

`true` iff `x ∈ L`.

### Algorithm

The point ``x`` belongs to the line if and only if ``a⋅x = b`` holds.
"""
function ∈(x::AbstractVector, L::Line2D)
    @assert length(x) == 2 "a $(length(x))-dimensional vector is " *
                           "incompatible with a 2-dimensional line"
    return _isapprox(dot(L.a, x), L.b)
end

"""
    rand(::Type{Line2D}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create a random line.

### Input

- `Line2D` -- type for dispatch
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
function rand(::Type{Line2D};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    @assert dim == 2 "cannot create a random Line2D of dimension $dim"
    rng = reseed!(rng, seed)
    a = randn(rng, N, dim)
    while iszero(a)
        a = randn(rng, N, dim)
    end
    b = randn(rng, N)
    return Line2D(a, b)
end

"""
    isempty(L::Line2D)

Check whether a 2D line is empty.

### Input

- `L` -- 2D line

### Output

`false`.
"""
function isempty(L::Line2D)
    return false
end

"""
    constrained_dimensions(L::Line2D)

Return the indices in which a 2D line is constrained.

### Input

- `L` -- 2D line

### Output

A vector of ascending indices `i` such that the line is constrained in dimension
`i`.

### Examples

A line with constraint ``x_i = 0`` (``i ∈ \\{1, 2\\}``) is only constrained in
dimension ``i``.
"""
function constrained_dimensions(L::Line2D)
    return nonzero_indices(L.a)
end

function _linear_map_hrep_helper(M::AbstractMatrix{N}, P::Line2D{N},
                                 algo::AbstractLinearMapAlgorithm) where {N}
    constraints = _linear_map_hrep(M, P, algo)
    if length(constraints) == 2
        # assuming these constraints define a line
        c = first(constraints)
        return Line2D(c.a, c.b)
    elseif isempty(constraints)
        return Universe{N}(size(M, 1))
    else
        error("unexpected number of $(length(constraints)) constraints")
    end
end

"""
    translate(L::Line2D, v::AbstractVector; [share]::Bool=false)

Translate (i.e., shift) a 2D line by a given vector.

### Input

- `L`     -- 2D line
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
function translate(L::Line2D, v::AbstractVector; share::Bool=false)
    @assert length(v) == dim(L) "cannot translate a $(dim(L))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    a = share ? L.a : copy(L.a)
    b = L.b + dot(L.a, v)
    return Line2D(a, b)
end

# the algorithm is a 2D specialization of the `Hyperplane` algorithm, except
# that it returns a `Singleton` for a 1D line
function project(L::Line2D{N}, block::AbstractVector{Int}; kwargs...) where {N}
    m = length(block)
    if m == 2
        @inbounds if block[1] == 1 && block[2] == 2
            return L  # no projection
        elseif block[1] == 2 && block[2] == 1
            return Line2D(L.a[block], L.b)  # swap a vector
        else
            throw(ArgumentError("invalid projection to $block"))
        end
    elseif m == 1
        # projection to dimension i
        cdims = constrained_dimensions(L)
        if length(cdims) == 1
            @inbounds if cdims[1] == block[1]
                # L: aᵢxᵢ = b where aᵢ ≠ 0
                return Singleton([L.b / L.a[cdims[1]]])
            else
                # L: aⱼxⱼ = b where i ≠ j
                return Universe{N}(1)
            end
        else
            # L is constrained in both dimensions
            @assert length(cdims) == 2
            return Universe{N}(1)
        end
    else
        throw(ArgumentError("cannot project a two-dimensional line to $m dimensions"))
    end
end

"""
    project(x::AbstractVector, L::Line2D)

Project a point onto a 2D line.

### Input

- `x` -- point/vector
- `L` -- 2D line

### Output

The projection of `x` onto `L`.

### Algorithm

The projection of ``x`` onto a line of the form ``a⋅x = b`` is

```math
    x - \\dfrac{a (a⋅x - b)}{‖a‖²}.
```
"""
function project(x::AbstractVector, L::Line2D)
    return x - L.a * (dot(L.a, x) - L.b) / norm(L.a, 2)^2
end
