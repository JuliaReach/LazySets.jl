import Base: rand,
             ∈,
             isempty

export Line,
       an_element,
       translate,
       translate!

"""
    Line{N, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}

Type that represents a line in the form

```math
    \\{y ∈ \\mathbb{R}^n: y = p + λn, λ ∈ \\mathbb{R}\\}
```
where ``p`` is a point on the line and ``n`` is its direction vector (not necessarily
normalized).

### Fields

- `p` -- point in the line
- `n` -- direction

### Examples

The line passing through the point ``[-1, 2, 3]`` and parallel to the vector
``[3, 0, -1]``:

```jldoctest
julia> Line([-1, 2, 3.], [3, 0, -1.])
Line{Float64,Array{Float64,1}}([-1.0, -1.0, 1.0], 1.0)
```
"""
struct Line{N, VN<:AbstractVector{N}} <: AbstractPolyhedron{N}
    p::VN
    n::VN

    # default constructor with length constraint
    function Line(p::VN, n::VN; check_direction::Bool=true) where {N, VN<:AbstractVector{N}}
        if check_direction
            @assert !iszero(n) "a line needs a non-zero direction vector"
        end
        return new{N, VN}(p, n)
    end
end

isoperationtype(::Type{<:Line}) = false
isconvextype(::Type{<:Line}) = true

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
direction(L::Line) = L.n

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
    normalize!(L.n, p)
    return Lx
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
    return Line(copy(L.p), normalize(L.n, p))
end

"""
    Line(p::AbstractVector, q::AbstractVector; [normalize]=true)

Constructor of a line give two points.

### Input

- `p`         -- point
- `q`         -- another point
- `normalize` -- (optional, default: `true`) if `true`, the direction of the line
                 has norm 1 (w.r.t the Euclidean norm)

### Output

The line which passes through `p` and `q`.

### Algorithm

Given two points ``p ∈ \\mathbb{R}^n`` and ``q ∈ \\mathbb{R}^n``, the line that
passes through these two points is
`L: `\\{y ∈ \\mathbb{R}^n: y = p + λ(q - p), λ ∈ \\mathbb{R}\\}``.
"""
function Line(p::AbstractVector, q::AbstractVector; normalize=true)
    d = q - p
    if normalize && iszero(d)
        throw(ArgumentError("points `p` and `q` should be distinct"))
    end
    n = d / dot(d, d)
    return Line(p, n)
end

# --- polyhedron interface functions ---

"""
    constraints_list(L::Line{N, VN}) where {N, VN}

Return the list of constraints of a line.

### Input

- `L` -- line

### Output

A list containing `n-1` half-spaces whose intersection is `L`, where `n` is the
ambient dimension of `L`.
"""
function constraints_list(L::Line{N, VN}) where {N, VN}
    p = L.p
    d = length(p)
    n = reshape(L.n, 1, d)
    K = nullspace(n)
    m = size(K, 2)
    @assert m == d - 1 "expected $(d - 1) normal half-spaces, got $m"

    out = Vector{HalfSpace{N, VN}}(undef, m)

    @inbounds for j in 1:m
        Kj = K[:, j]
        b = dot(Kj, p)
        out[j] = HalfSpace(Kj, b)
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
function ρ(d::AbstractVector, L::Line)
    if isapproxzero(d, L.n)
        return dot(d, L.p)
    else
        return Inf
    end
end

"""
    σ(d::AbstractVector{N}, L::Line{N}) where {N<:Real}

Return the support vector of a line in a given direction.

### Input

- `d` -- direction
- `L` -- line

### Output

The support vector in the given direction.
"""
function σ(d::AbstractVector, L::Line)
    if isapproxzero(d, L.n)
        return L.p
    else
        error("the support vector is undefined because the line is unbounded")
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

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ P``
"""
isuniversal(L::Line; witness::Bool=false) = isuniversal(L, Val(witness))

isuniversal(L::Line, ::Val{false}) = false

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

The point ``x`` belongs to the line ``L : p + λ⋅n`` if and only if
``x - p`` is proportional to the direction ``n``.
"""
function ∈(x::AbstractVector, L::Line)
    @assert length(x) == dim(L) "expected the point and the line to have the same dimension, " *
                                "but they are $(length(x)) and $(dim(L)) respectively"
    return samedir(x - L.p, L.n)
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
    n = randn(rng, N, dim)
    while iszero(n)
        n = randn(rng, N, dim)
    end
    p = randn(rng, N, dim)
    return Line(p, n)
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
    constrained_dimensions(L::Line)

Return the indices in which a line is constrained.

### Input

- `L` -- line

### Output

A vector of ascending indices `i` such that the line is constrained in dimension
`i`.

### Algorithm

For each coordinate ``i``, vector of the form ``x_i = p_i + λ n_i`` are constrained
i.e. (they belong to the line and are bounded) if and only if ``n_i`` is zero.
Hence, this function returns all indices of the normal vector ``n`` for which
the ``i``-th coordinate is nonzero.

### Examples

The line ``y = 5`` in two dimensions can be written as ``p = [0, 5]`` and
``n = [1, 0)]``. This line constrains dimension ``2``.

```jldoctest
julia> constrained_dimensions(Line([0, 5.], [1, 0.]))
1-element Array{Int64,1}:
 2
```
"""
function constrained_dimensions(L::Line)
    return findall(isapproxzero, L.n)
end

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
    @assert length(v) == dim(L) "cannot translate a $(dim(L))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
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
    L.p .+= v
    return L
end
