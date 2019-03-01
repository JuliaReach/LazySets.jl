"""
    decompose(S::LazySet{N},
              partition::AbstractVector{<:AbstractVector{Int}},
              block2oa
             )::CartesianProductArray{N} where {N<:Real}

Decompose a high-dimensional set into a Cartesian product of overapproximations
of the projections over the specified subspaces.

### Input

- `S`         -- set
- `partition` -- vector of blocks (i.e., of vectors of integers) (see the Notes
                 below)
- `block2oa`  -- mapping from block indices in `partition` to a corresponding
                 overapproximation option; we only require access via `[⋅]` (but
                 see also the Notes below)

### Output

A `CartesianProductArray` containing the low-dimensional approximated
projections.

### Algorithm

For each block a specific `project` method is called, dispatching on the
corresponding overapproximation option.

### Notes

The argument `partition` requires some discussion.
Typically, the list of blocks should form a partition of the set
``\\{1, \\dots, n\\}`` represented as a list of consecutive blocks, where ``n``
is the ambient dimension of set `S`.

However, technically there is no problem if the blocks are not consecutive,
blocks are missing, blocks occur more than once, or blocks are overlapping.
This function will, however, stick to the order of blocks, so the resulting set
must be interpreted with care in such cases.
One use case is the need of a projection consisting of several blocks.

For convenience, the argument `block2oa` can also be given as a single `oa`
option, which is then interpreted as the option for all blocks.
"""
function decompose(S::LazySet{N},
                   partition::AbstractVector{<:AbstractVector{Int}},
                   block2oa
                  )::CartesianProductArray{N} where {N<:Real}
    n = dim(S)
    result = Vector{LazySet{N}}(undef, length(partition))

    @inbounds for (i, block) in enumerate(partition)
        result[i] = project(S, block, block2oa[i], n)
    end
    return CartesianProductArray(result)
end

# convenience method
function decompose(S::LazySet{N},
                   partition::AbstractVector{<:AbstractVector{Int}},
                   oa::Union{Type{<:LazySet},
                             Pair{<:UnionAll, <:Real},
                             Real,
                             Type{<:AbstractDirections}
                            }
                  )::CartesianProductArray{N} where {N<:Real}
    n = dim(S)
    result = Vector{LazySet{N}}(undef, length(partition))

    @inbounds for (i, block) in enumerate(partition)
        result[i] = project(S, block, oa, n)
    end
    return CartesianProductArray(result)
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            set_type::Type{<:LinearMap},
            [n]::Int=dim(S)
           )::LinearMap{N} where {N<:Real}

Project a high-dimensional set to a given block by using a lazy linear map.

### Input

- `S`         -- set
- `block`     -- block structure - a vector with the dimensions of interest
- `LinearMap` -- used for dispatch
- `n`         -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A lazy `LinearMap` representing a projection of the set `S` to block `block`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         set_type::Type{<:LinearMap},
                         n::Int=dim(S)
                        )::LinearMap{N} where {N<:Real}
    m = length(block)
    M = sparse(1:m, block, ones(N, m), m, n)
    return M * S
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            set_type::Type{<:LazySet},
            [n]::Int=dim(S)
           ) where {N<:Real}

Project a high-dimensional set to a given block and set type, possibly involving
an overapproximation.

### Input

- `S`        -- set
- `block`    -- block structure - a vector with the dimensions of interest
- `set_type` -- target set type
- `n`        -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A set of type `set_type` representing an overapproximation of the projection of
`S`.

### Algorithm

1. Project the set `S` with `M⋅S`, where `M` is the identity matrix in the block
coordinates and zero otherwise.
2. Overapproximate the projected lazy set using `overapproximate` and
`set_type`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         set_type::Type{<:LazySet},
                         n::Int=dim(S)
                        ) where {N<:Real}
    lm = project(S, block, LinearMap, n)
    return overapproximate(lm, set_type)
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            set_type_and_precision::Pair{<:UnionAll, <:Real},
            [n]::Int=dim(S)
           ) where {N<:Real}

Project a high-dimensional set to a given block and set type with a certified
error bound.

### Input

- `S`     -- set
- `block` -- block structure - a vector with the dimensions of interest
- `set_type_and_precision` -- pair `(T, ε)` of a target set type `T` and an
                              error bound `ε` for approximation
- `n`     -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A set representing the epsilon-close approximation of the projection of `S`.

### Notes

Currently we only support `HPolygon` as set type, which implies that the set
must be two-dimensional.

### Algorithm

1. Project the set `S` with `M⋅S`, where `M` is the identity matrix in the block
coordinates and zero otherwise.
2. Overapproximate the projected lazy set with the given error bound `ε`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         set_type_and_precision::Pair{<:UnionAll, <:Real},
                         n::Int=dim(S)
                        ) where {N<:Real}
    set_type = set_type_and_precision[1]
    ε = set_type_and_precision[2]
    @assert length(block) == 2 && set_type == HPolygon "currently only 2D " *
        "HPolygon decomposition is supported"

    lm = project(S, block, LinearMap, n)
    return overapproximate(lm, set_type, ε)
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            ε::Real,
            [n]::Int=dim(S)
           ) where {N<:Real}

Project a high-dimensional set to a given block and set type with a certified
error bound.

### Input

- `S`     -- set
- `block` -- block structure - a vector with the dimensions of interest
- `ε`     -- error bound for approximation
- `n`     -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

A set representing the epsilon-close approximation of the projection of `S`.

### Algorithm

1. Project the set `S` with `M⋅S`, where `M` is the identity matrix in the block
coordinates and zero otherwise.
2. Overapproximate the projected lazy set with the given error bound `ε`.
The target set type is chosen automatically.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         ε::Real,
                         n::Int=dim(S)
                        ) where {N<:Real}
    # currently we only support HPolygon
    if length(block) == 2
        set_type = HPolygon
    else
        throw(ArgumentError("ε-close approximation is only supported for 2D " *
                            "blocks"))
    end
    return project(S, block, set_type => ε, n)
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            directions::Type{<:AbstractDirections},
            [n]::Int
           ) where {N<:Real}

Project a high-dimensional set to a given block using template directions.

### Input

- `S`          -- set
- `block`      -- block structure - a vector with the dimensions of interest
- `directions` -- template directions
- `n`          -- (optional, default: `dim(S)`) ambient dimension of the set `S`

### Output

The template direction approximation of the projection of `S`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         directions::Type{<:AbstractDirections},
                         n::Int=dim(S)
                        ) where {N<:Real}
    lm = project(S, block, LinearMap, n)
    return overapproximate(lm, directions(length(block)))
end
