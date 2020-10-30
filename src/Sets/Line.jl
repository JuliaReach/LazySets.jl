import Base: rand,
             ∈,
             isempty

export Line,
       an_element,
       translate,
       translate!,
       distance

"""
    Line{N, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a line in the form

```math
    \\{y ∈ \\mathbb{R}^n: y = p + λd, λ ∈ \\mathbb{R}\\}
```
where ``p`` is a point on the line and ``d`` is its direction vector (not necessarily
normalized).

### Fields

- `p` -- point on the line
- `d` -- direction

### Examples

The line passing through the point ``[-1, 2, 3]`` and parallel to the vector
``[3, 0, -1]``:

```jldoctest
julia> Line([-1, 2, 3.], [3, 0, -1.])
Line{Float64,Array{Float64,1}}([-1.0, 2.0, 3.0], [3.0, 0.0, -1.0])
```
"""
struct Line{N, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    p::VN
    d::VN

    # default constructor with length constraint
    function Line(p::VN, d::VN; check_direction::Bool=true) where {N, VN<:AbstractVector{N}}
        if check_direction
            @assert !iszero(d) "a line needs a non-zero direction vector"
        end
        return new{N, VN}(p, d)
    end
end

isoperationtype(::Type{<:Line}) = false
isconvextype(::Type{<:Line}) = true

"""
    Line(; from::AbstractVector, to::AbstractVector, normalize=true)

Constructor of a line given two points.

### Input

- `from`      -- point
- `to`        -- another point
- `normalize` -- (optional, default: `true`) if `true`, the direction of the line
                 has norm 1 (w.r.t the Euclidean norm)

### Output

The line which passes through `from` and `to`.

### Algorithm

Given two points ``p ∈ \\mathbb{R}^n`` and ``q ∈ \\mathbb{R}^n``, the line that
passes through these two points is
`L: `\\{y ∈ \\mathbb{R}^n: y = p + λ(q - p), λ ∈ \\mathbb{R}\\}``.
"""
function Line(; from::AbstractVector, to::AbstractVector, normalize=true)
    d = from - to
    if normalize && iszero(d)
        throw(ArgumentError("points `$from` and `$to` should be distinct"))
    end
    d_n = normalize ? d / dot(d, d) : d
    return Line(from, d_n)
end

"""
    Line(a::AbstractVector{N}, b::N; normalize=true) where {N}

Constructor of a two-dimensional line given `ax = b`.

### Input

- `a`         -- two-dimensional vector
- `b`         -- scalar
- `normalize` -- (optional, default: `true`) if `true`, the direction of the line
                 has norm 1 (w.r.t the Euclidean norm)

### Output

The `Line` that satisfies `a ⋅ x = b`.
"""
function Line(a::AbstractVector{N}, b::N; normalize=true) where {N}
    @assert length(a) == 2 "expected a normal vector of length two, but it is $(length(a))-dimensional"

    got_horizontal = iszero(a[1])
    got_vertical = iszero(a[2])

    if got_horizontal && got_vertical
        throw(ArgumentError("the vector $a must be non-zero"))
    end

    if got_horizontal
        α = b / a[2]
        p = [zero(N), α]
        q = [one(N), α]
    elseif got_vertical
        β = b / a[1]
        p = [β, zero(N)]
        q = [β, one(N)]
    else
        α = b / a[2]
        μ = a[1] / a[2]
        p = [zero(N), α]
        q = [one(N), α - μ]
    end
    return Line(from=p, to=q, normalize=normalize)
end

"""
    direction(L::Line)

Return the direction of the line.

### Input

- `L` -- line

### Output

The line's field corresponding to the direction to the line.

### Notes

The direction is not necessarily normalized.
See [`normalize(::Line, ::Real)`](@ref) / [`normalize!(::Line, ::Real)`](@ref)
for such operation.
"""
direction(L::Line) = L.d

"""
    normalize!(L::Line, p::Real=2.0)

Normalize the direction of a line storing the result in `L`.

### Input

- `L` -- line
- `p` -- (optional, default: `2.0`) vector `p`-norm used in the normalization

### Output

A line whose direction has unit norm w.r.t the given `p`-norm.
"""
function normalize!(L::Line, p::Real=2.0)
    normalize!(L.d, p)
    return L
end

"""
    normalize(L::Line, p::Real=2.0)

Normalize the direction of a line.

### Input

- `L` -- line
- `p` -- (optional, default: `2.0`) vector `p`-norm used in the normalization

### Output

A line whose direction has unit norm w.r.t the given `p`-norm.

### Notes

See also [`normalize!(::Line, ::Real)`](@ref) for the in-place version.
"""
function normalize(L::Line, p::Real=2.0)
    normalize!(copy(L), p)
end

# --- polyhedron interface functions ---

"""
    constraints_list(L::Line{N, VN}) where {N, VN}

Return the list of constraints of a line.

### Input

- `L` -- line

### Output

A list containing `2n-2` half-spaces whose intersection is `L`, where `n` is the
ambient dimension of `L`.
"""
function constraints_list(L::Line{N, VN}) where {N, VN}
    p = L.p
    n = length(p)
    d = reshape(L.d, 1, n)
    K = nullspace(d)
    m = size(K, 2)
    @assert m == n - 1 "expected $(n - 1) normal half-spaces, got $m"

    out = Vector{HalfSpace{N, VN}}(undef, 2m)
    idx = 1
    @inbounds for j in 1:m
        Kj = K[:, j]
        b = dot(Kj, p)
        out[idx] = HalfSpace(Kj, b)
        out[idx+1] = HalfSpace(-Kj, -b)
        idx += 2
    end
    return out
end

# --- LazySet interface functions ---


"""
    dim(L::Line)

Return the ambient dimension of a line.

### Input

- `L` -- line

### Output

The ambient dimension of the line.
"""
dim(L::Line) = length(L.p)

"""
    ρ(d::AbstractVector, L::Line)

Return the support function of a line in a given direction.

### Input

- `d` -- direction
- `L` -- line

### Output

The support function in the given direction.
"""
ρ(d::AbstractVector, L::Line) = _ρ(d, L)
ρ(d::AbstractVector{N}, L::Line{N, <:AbstractVector{N}}) where {N<:Real} = _ρ(d, L) # disambiguation

function _ρ(d::AbstractVector, L::Line)
    if isapproxzero(dot(d, L.d))
        return dot(d, L.p)
    else
        return Inf
    end
end

"""
    σ(d::AbstractVector, L::Line)

Return the support vector of a line in a given direction.

### Input

- `d` -- direction
- `L` -- line

### Output

The support vector in the given direction.
"""
function σ(d::AbstractVector, L::Line)
    if isapproxzero(dot(d, L.d))
        return L.p
    else
        throw(ArgumentError("the support vector is undefined because the line is unbounded"))
    end
end

"""
    isbounded(L::Line)

Determine whether a line is bounded.

### Input

- `L` -- line

### Output

`false`.
"""
isbounded(::Line) = false

"""
    isuniversal(L::Line; [witness::Bool]=false)

Check whether a line is universal.

### Input

- `P`       -- line
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` is `false`: `true` if the ambient dimension is one, `false` otherwise

* If `witness` is `true`: `(true, [])` if the ambient dimension is one,
  `(false, v)` where ``v ∉ P`` otherwise
"""
isuniversal(L::Line; witness::Bool=false) = isuniversal(L, Val(witness))

# TODO implement case with witness
isuniversal(L::Line, ::Val{false}) = dim(L) == 1

"""
    an_element(L::Line)

Return some element of a line.

### Input

- `L` -- line

### Output

An element on the line.
"""
an_element(L::Line) = L.p

"""
    ∈(x::AbstractVector, L::Line)

Check whether a given point is contained in a line.

### Input

- `x` -- point/vector
- `L` -- line

### Output

`true` iff `x ∈ L`.

### Algorithm

The point ``x`` belongs to the line ``L : p + λd`` if and only if
``x - p`` is proportional to the direction ``d``.
"""
∈(x::AbstractVector, L::Line) = __in(x, L)
∈(x::AbstractVector{N}, L::Line{N, VN}) where {N<:Real, VN<:AbstractVector{N}} = __in(x, L)

function __in(x::AbstractVector, L::Line)
    @assert length(x) == dim(L) "expected the point and the line to have the same dimension, " *
                                "but they are $(length(x)) and $(dim(L)) respectively"
    _isapprox(x, L.p) && return true

    # TODO pass "minus" option to samedir
    return first(samedir(x - L.p, L.d)) || first(samedir(x - L.p, -L.d))
end

"""
    rand(::Type{Line}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

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
"""
function rand(::Type{Line};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing)

    rng = reseed(rng, seed)
    d = randn(rng, N, dim)
    while iszero(d)
        d = randn(rng, N, dim)
    end
    p = randn(rng, N, dim)
    return Line(p, d)
end

"""
    isempty(L::Line)

Return if a line is empty or not.

### Input

- `L` -- line

### Output

`false`.
"""
isempty(::Line) = false

"""
    translate(L::Line, v::AbstractVector)

Translate (i.e., shift) a line by a given vector.

### Input

- `L` -- line
- `v` -- translation vector

### Output

A translated line.

### Notes

See also `translate!` for the in-place version.
"""
function translate(L::Line, v::AbstractVector)
    return translate!(copy(L), v)
end

"""
    translate!(L::Line, v::AbstractVector)

Translate (i.e., shift) a line by a given vector storing the result in `L`.

### Input

- `L` -- line
- `v` -- translation vector

### Output

A translated line, modifying `L` in-place.
"""
function translate!(L::Line, v::AbstractVector)
    @assert length(v) == dim(L) "cannot translate a $(dim(L))-dimensional " *
                                "set by a $(length(v))-dimensional vector"

    L.p .+= v
    return L
end

"""
    distance(x::AbstractVector, L::Line, p::Real=2.0)

Compute the distance between point `x` and the line with respect to the given
`p`-norm.

### Input

- `x` -- vector
- `L` -- line
- `p` -- (optional, default: `2.0`) the `p`-norm used; `p = 2.0` corresponds to
         the usual Euclidean norm

### Output

A scalar representing the distance between `x` and the line `L`.
"""
function distance(x::AbstractVector, L::Line, p::Real=2.0)
    d = L.d  # direction of the line
    t = dot(x - L.p, d) / dot(d, d)
    return distance(x, L.p + t*d, p)
end
distance(L::Line, x::AbstractVector, p::Real=2.0) = distance(x, L, p)

"""
    distance(x::AbstractSingleton, L::Line, p::Real=2.0)

Compute the distance between the singleton `x` and the line with respect to the given
`p`-norm.

### Input

- `x` -- singleton, i.e. a set with one element
- `L` -- line
- `p` -- (optional, default: `2.0`) the `p`-norm used; `p = 2.0` corresponds to
         the usual Euclidean norm

### Output

A scalar representing the distance between the element wrapped by `x` and the line `L`.
"""
function distance(x::AbstractSingleton, L::Line, p::Real=2.0)
    return distance(element(x), L, p)
end
distance(L::Line, x::AbstractSingleton, p::Real=2.0) = distance(x, L, p)

"""
    linear_map(M::AbstractMatrix, L::Line)

Concrete linear map of a line.

### Input

- `M` -- matrix
- `L` -- line

### Output

The line obtained by applying the linear map to the point and direction of `L`.
"""
linear_map(M::AbstractMatrix, L::Line) = _linear_map(M, L)

linear_map(M::AbstractMatrix{N}, L::Line{N,VN}) where {N<:Real, VN<:AbstractVector{N}} = _linear_map(M, L)

function _linear_map(M::AbstractMatrix, L::Line)
    @assert dim(L) == size(M, 2) "a linear map of size $(size(M)) cannot be " *
                                 "applied to a set of dimension $(dim(L))"

    Mp = M * L.p
    Md = M * L.d
    return Line(Mp, Md)
end
