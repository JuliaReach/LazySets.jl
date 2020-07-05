import Base: rand,
             ∈,
             isempty

export Line, an_element

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
    function Line(p::VN, n::N; check_direction::Bool=true) where {N, VN<:AbstractVector{N}}
        if check_direction
            @assert !iszero(n) "a line needs a non-zero direction vector"
        end
        return new{N, VN}(p, n)
    end
end

isoperationtype(::Type{<:Line}) = false
isconvextype(::Type{<:Line}) = true

direction(L::Line) = L.n
normalize!(L::Line, p::Real=2) = normalize!(L.n, p)
normalize(L::Line, p::Real=2) = normalize!(copy(L), p)

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
`L: `\\{y ∈ \\mathbb{R}^n: y = p + λ(p - q), λ ∈ \\mathbb{R}\\}``.
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
    return _constraints_list_hyperplane(L.a, L.b)
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

### Algorithm

Witness production falls back to `isuniversal(::Hyperplane)`.
"""
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
    p = randn(rng, N)
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

### Examples

A line with constraint ``x1 = 0`` is constrained in dimension 1 only.
"""
function constrained_dimensions(L::Line)
    # TODO
end


function _linear_map_hrep_helper(M::AbstractMatrix{N}, P::Line{N},
                                 algo::AbstractLinearMapAlgorithm) where {N<:Real}
    # TODO
    #=
    constraints = _linear_map_hrep(M, P, algo)
    if length(constraints) == 2
        # assuming these constraints define a line
        c = first(constraints)
        return Line(c.a, c.b)
    elseif isempty(constraints)
        return Universe{N}(size(M, 1))
    else
        error("unexpected number of $(length(constraints)) constraints")
    end
    =#
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
    return translate(copy(L), v)
end

"""
    translate!(L::Line, v::AbstractVector)

Translate (i.e., shift) a line by a given vector in-place.

### Input

- `L` -- line
- `v` -- translation vector

### Output

A translated line.
"""
function translate!(L::Line, v::AbstractVector)
    L.p .+= v
    return L
end
