export EmptySet, ∅

"""
    EmptySet{N} <: ConvexSet{N}

Type that represents the empty set, i.e., the set with no elements.
"""
struct EmptySet{N} <: ConvexSet{N}
    dim::Int
end

# default constructor of type Float64
EmptySet(n::Int) = EmptySet{Float64}(n)

isoperationtype(::Type{<:EmptySet}) = false
isconvextype(::Type{<:EmptySet}) = true

"""
    ∅

Alias for `EmptySet{Float64}`.
"""
const ∅ = EmptySet{Float64}

"""
    dim(∅::EmptySet)

Return the dimension of an empty set.

### Input

- `∅` -- an empty set

### Output

The dimension of the empty set.
"""
function dim(∅::EmptySet)
    return ∅.dim
end

"""
    σ(d::AbstractVector, ∅::EmptySet)

Return the support vector of an empty set.

### Input

- `d` -- direction
- `∅` -- empty set

### Output

An error.
"""
function σ(d::AbstractVector, ∅::EmptySet)
    return error("the support vector of an empty set is undefined")
end

"""
    ρ(d::AbstractVector, ∅::EmptySet)

Evaluate the support function of an empty set in a given direction.

### Input

- `d` -- direction
- `∅` -- empty set

### Output

An error.
"""
function ρ(d::AbstractVector, ∅::EmptySet)
    return error("the support function of an empty set is undefined")
end

function isboundedtype(::Type{<:EmptySet})
    return true
end

"""
    isbounded(∅::EmptySet)

Check whether an empty set is bounded.

### Input

- `∅` -- empty set

### Output

`true`.
"""
function isbounded(::EmptySet)
    return true
end

"""
    isuniversal(∅::EmptySet{N}, [witness]::Bool=false) where {N}

Check whether an empty set is universal.

### Input

- `∅`       -- empty set
- `witness` -- (optional, default: `false`) compute a witness if activated

### Output

* If `witness` option is deactivated: `false`
* If `witness` option is activated: `(false, v)` where ``v ∉ S``
"""
function isuniversal(∅::EmptySet{N}, witness::Bool=false) where {N}
    if witness
        return (false, zeros(N, dim(∅)))
    else
        return false
    end
end

"""
    ∈(x::AbstractVector, ∅::EmptySet)

Check whether a given point is contained in an empty set.

### Input

- `x` -- point/vector
- `∅` -- empty set

### Output

`false`.

### Examples

```jldoctest
julia> [1.0, 0.0] ∈ ∅(2)
false
```
"""
function ∈(x::AbstractVector, ∅::EmptySet)
    return false
end

"""
    an_element(∅::EmptySet)

Return some element of an empty set.

### Input

- `∅` -- empty set

### Output

An error.
"""
function an_element(∅::EmptySet)
    return error("an empty set does not contain any element")
end

"""
    rand(::Type{EmptySet}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

Create an empty set (note that there is nothing to randomize).

### Input

- `EmptySet` -- type for dispatch
- `N`        -- (optional, default: `Float64`) numeric type
- `dim`      -- (optional, default: 2) dimension
- `rng`      -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`     -- (optional, default: `nothing`) seed for reseeding

### Output

The (only) empty set of the given numeric type and dimension.
"""
function rand(::Type{EmptySet};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    return EmptySet{N}(dim)
end

"""
    isempty(∅::EmptySet)

Check if the empty set is empty.

### Input

- `∅` -- empty set

### Output

`true`.
"""
function isempty(∅::EmptySet)
    return true
end

"""
    norm(∅::EmptySet, [p]::Real=Inf)

Return the norm of an empty set.
It is the norm of the enclosing ball (of the given ``p``-norm) of minimal volume
that is centered in the origin.

### Input

- `∅` -- empty set
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function norm(∅::EmptySet, p::Real=Inf)
    return error("the norm of an empty set is undefined")
end

"""
    radius(∅::EmptySet, [p]::Real=Inf)

Return the radius of an empty set.
It is the radius of the enclosing ball (of the given ``p``-norm) of minimal
volume with the same center.

### Input

- `∅` -- empty set
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function radius(∅::EmptySet, p::Real=Inf)
    return error("the radius of an empty set is undefined")
end

"""
    diameter(∅::EmptySet, [p]::Real=Inf)

Return the diameter of an empty set.
It is the maximum distance between any two elements of the set or, equivalently,
the diameter of the enclosing ball (of the given ``p``-norm) of minimal volume
with the same center.

### Input

- `∅` -- empty set
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function diameter(∅::EmptySet, p::Real=Inf)
    return error("the diameter of an empty set is undefined")
end

"""
    vertices(∅::EmptySet)

Construct an iterator over the vertices of an empty set.

### Input

- `∅` -- empty set

### Output

The empty iterator, as the empty set does not contain any vertices.
"""
function vertices(∅::EmptySet)
    N = eltype(∅)
    return EmptyIterator{Vector{N}}()
end

"""
    vertices_list(∅::EmptySet)

Return the list of vertices of an empty set.

### Input

- `∅` -- empty set

### Output

The empty list of vertices, as the empty set does not contain any vertices.
"""
function vertices_list(∅::EmptySet)
    N = eltype(∅)
    return Vector{Vector{N}}()
end

"""
    linear_map(M::AbstractMatrix{N}, ∅::EmptySet{N}) where {N}

Return the linear map of an empty set.

### Input

- `M` -- matrix
- `∅` -- empty set

### Output

An empty set.
"""
function linear_map(M::AbstractMatrix, ∅::EmptySet)
    N = eltype(∅)
    @assert size(M, 2) == dim(∅) "cannot apply a $(size(M))-dimensional " *
                                 "matrix to a $(dim(∅))-dimensional set"

    return EmptySet{N}(size(M, 1))
end

"""
    translate(∅::EmptySet, v::AbstractVector)

Translate (i.e., shift) an empty set by a given vector.

### Input

- `∅` -- empty set
- `v` -- translation vector

### Output

The empty set.
"""
function translate(∅::EmptySet, v::AbstractVector)
    return translate!(∅, v)  # no need to copy
end

function translate!(∅::EmptySet, v::AbstractVector)
    @assert length(v) == dim(∅) "cannot translate a $(dim(∅))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return ∅
end

"""
    plot_recipe(∅::EmptySet{N}, [ε]=zero(N)) where {N}

Convert an empty set to a sequence of points for plotting.
In the special case of an empty set, the sequence is empty.

### Input

- `∅` -- empty set
- `ε` -- (optional, default: `0`) ignored, used for dispatch

### Output

An empty array.
"""
function plot_recipe(∅::EmptySet{N}, ε=zero(N)) where {N}
    return []
end

"""
    area(∅::EmptySet)

Return the area of an empty set.

### Input

- `∅` -- empty set

### Output

``0``.
"""
function area(∅::EmptySet)
    N = eltype(∅)
    return zero(N)
end

"""
    volume(∅::EmptySet{N}) where {N}

Return the volume of an empty set.

### Input

- `∅` -- empty set

### Output

``0``.
"""
function volume(∅::EmptySet{N}) where {N}
    return zero(N)
end

function project(∅::EmptySet{N}, block::AbstractVector{Int}; kwargs...) where {N}
    return EmptySet{N}(length(block))
end

"""
    chebyshev_center_radius(∅::EmptySet; [kwargs]...)

Compute the [Chebyshev center](https://en.wikipedia.org/wiki/Chebyshev_center)
and the corresponding radius of an empty set.

### Input

- `∅`      -- empty set
- `kwargs` -- further keyword arguments (ignored)

### Output

An error.
"""
function chebyshev_center_radius(∅::EmptySet; kwargs...)
    return error("the Chebyshev center and radius of an empty set are undefined")
end

"""
    low(∅::EmptySet)

Return a vector with the lowest coordinates of an empty set in each canonical
direction.

### Input

- `∅` -- empty set

### Output

An error.

### Notes

See also [`low(∅::EmptySet, i::Int)`](@ref).
"""
function low(∅::EmptySet)
    return error("the lower bound of an empty set is undefined")
end

"""
    high(∅::EmptySet)

Return a vector with the highest coordinates of an empty set in each canonical
direction.

### Input

- `∅` -- empty set

### Output

An error.

### Notes

See also [`high(∅::EmptySet, i::Int)`](@ref).
"""
function high(∅::EmptySet)
    return error("the upper bound of an empty set is undefined")
end

"""
    low(∅::EmptySet, i::Int)

Return the lowest coordinate of an empty set in the given direction.

### Input

- `∅` -- empty set
- `i` -- dimension of interest

### Output

An error.
"""
function low(∅::EmptySet, i::Int)
    return error("the lower bound of an empty set is undefined")
end

"""
    high(∅::EmptySet, i::Int)

Return the highest coordinate of an empty set in the given direction.

### Input

- `∅` -- empty set
- `i` -- dimension of interest

### Output

An error.
"""
function high(∅::EmptySet, i::Int)
    return error("the upper bound of an empty set is undefined")
end

"""
    rectify(∅::EmptySet)

Concrete rectification of an empty set.

### Input

- `∅` -- empty set

### Output

The empty set.
"""
function rectify(∅::EmptySet)
    return ∅
end

"""
    complement(∅::EmptySet{N}) where {N}

Return the complement of an empty set.

### Input

- `∅` -- empty set

### Output

The universe of the same dimension.
"""
function complement(∅::EmptySet{N}) where {N}
    return Universe{N}(dim(∅))
end

"""
    reflect(∅::EmptySet)

Concrete reflection of an empty set.

### Input

- `∅` -- empty set

### Output

The same empty set.
"""
function reflect(∅::EmptySet)
    return ∅
end

function scale(::Real, ∅::EmptySet)
    return ∅
end

function scale!(::Real, ∅::EmptySet)
    return ∅
end
