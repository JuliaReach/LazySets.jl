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

direction(L::Line) = L.n

function normalize!(L::Line, p::Real=2)
    normalize!(L.n, p)
    return L
end

function normalize(L::Line, p::Real=2.)
    return Line(copy(L.p), normalize(L.n, p))
end

function distance(x::AbstractVector, y::AbstractVector, p::Real=2.)
    return norm(x - y, p)
end

function distance(x::AbstractVector, L::Line, p::Real=2.)
    q = L.p  # point in the line
    n = L.n  # direction of the line

    t = dot(x - q, n) / dot(n, n)
    return distance(x, q + t*n, p)
end

"""
    Line(p::AbstractVector, q::AbstractVector)

Constructor of a line give two points.

### Input

- `p` -- point
- `q` -- another point

### Output

The line which passes through `p` and `q`.

### Algorithm

Given two points ``p ∈ \\mathbb{R}^n`` and ``q ∈ \\mathbb{R}^n``, the line that
passes through these two points is
`L: `\\{y ∈ \\mathbb{R}^n: y = p + λ(q - p), λ ∈ \\mathbb{R}\\}``.
"""
function Line(p::AbstractVector, q::AbstractVector)
    Line(p, q - p)
end

# --- polyhedron interface functions ---

"""
    constraints_list(L::Line)

Return the list of constraints of a line.

### Input

- `L` -- line

### Output

A list containing two half-spaces.
"""
function constraints_list(L::Line)
    # TODO
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
        return Inf
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
