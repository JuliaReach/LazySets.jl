"""
    default_block_structure(S::LazySet, set_type::Type{<:LazySet})::AbstractVector{Int}

Compute the default block structure.

### Input

- `S`        -- set
- `set_type` -- target set type

### Output

A vector representing the block structure, such that:

- If the target `set_type` is an interval, the default is blocks of size 1.
- Otherwise, the default is blocks of size 2. Depending on the dimension,
  the last block has size 1 or 2.
"""
@inline function default_block_structure(S::LazySet, set_type::Type{<:LazySet})::AbstractVector{Int}
    n = dim(S)
    if set_type == Interval
        return fill(1, n)
    else
        if n % 2 == 0
            return fill(2, div(n, 2))
        else
            res = fill(2, div(n+1, 2))
            res[end] = 1
            return res
        end
    end
end

"""
    decompose(S::LazySet{N};
              [set_type]::Type{<:Union{HPolygon, Hyperrectangle, LazySets.Interval}}=Hyperrectangle,
              [ɛ]::Real=Inf,
              [blocks]::AbstractVector{Int}=default_block_structure(S, set_type),
             )::CartesianProductArray where {N<:Real}

Decompose a high-dimensional set into a Cartesian product of overapproximations
of the projections over the specified subspaces.

### Input

- `S`        -- set
- `set_type` -- (optional, default: `Hyperrectangle`) type of set approximation
                for each subspace
- `ɛ`        -- (optional, default: `Inf`) error bound for polytopic approximation
- `blocks`   -- (optional, default: [2, …, 2] or [1, …, 1] if `set_type` is an interval)
                block structure - a vector with the size of each block

### Output

A `CartesianProductArray` containing the low-dimensional approximated
projections.

### Algorithm

For each block a specific `project` method is called, dispatched on the
`set_type` argument.
"""
function decompose(S::LazySet{N};
                   set_type::Type{<:Union{HPolygon, Hyperrectangle, LazySets.Interval}}=Hyperrectangle,
                   ɛ::Real=Inf,
                   blocks::AbstractVector{Int}=default_block_structure(S, set_type)
                  )::CartesianProductArray where {N<:Real}
    n = dim(S)
    result = Vector{set_type{N}}()
    block_start = 1
    @inbounds for bi in blocks
        push!(result,
              project(S, block_start:(block_start + bi - 1), set_type, n, ɛ))
        block_start += bi
    end
    return CartesianProductArray(result)
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            set_type::Type{<:LazySet},
            [n]::Int=dim(S),
            [ɛ]::Real=Inf
           )::LazySet{N} where {N<:Real}

Default implementation for projecting a high-dimensional set to a given set type
with possible overapproximation.

### Input

- `S` -- set
- `block` -- block structure - a vector with the dimensions of interest
- `set_type` -- target set type
- `n` -- (optional, default: `dim(S)`) ambient dimension of the set `S`
- `ɛ` -- (optional, default: `Inf`) ignored

### Output

A set of type `set_type` representing an overapproximation of the projection of
`S`.

### Algorithm

1. Project the set `S` with `M⋅S`, where `M` is the identity matrix in the block coordinates and zero otherwise.
2. Overapproximate the projected lazy set using `overapproximate`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         set_type::Type{<:LazySet},
                         n::Int=dim(S),
                         ɛ::Real=Inf
                        )::LazySet{N} where {N<:Real}
    M = sparse(1:length(block), block, ones(N, length(block)), length(block), n)
    return overapproximate(M * S, set_type)
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            set_type::Type{<:HPolygon},
            [n]::Int=dim(S),
            [ɛ]::Real=Inf
           )::HPolygon where {N<:Real}

Project a high-dimensional set to a two-dimensional polygon with a certified
error bound.

### Input

- `S` -- set
- `block` -- block structure - a vector with the dimensions of interest
- `set_type` -- `HPolygon` - used for dispatch
- `n` -- (optional, default: `dim(S)`) ambient dimension of the set `S`
- `ɛ` -- (optional, default: `Inf`) error bound for polytopic approximation

### Output

A `HPolygon` representing the epsilon-close approximation of the box
approximation of the projection of `S`.

### Notes

`block` must have length 2.

### Algorithm

If `ɛ < Inf`, the algorithm proceeds as follows:

1. Project the set `S` with `M⋅S`, where `M` is the identity matrix in the block coordinates and zero otherwise.
2. Overapproximate the set with the given error bound `ɛ`.

If `ɛ == Inf`, the algorithm uses a box approximation.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         set_type::Type{<:HPolygon},
                         n::Int=dim(S),
                         ɛ::Real=Inf
                        )::HPolygon where {N<:Real}
    @assert length(block) == 2 "only 2D HPolygon decomposition is supported"

    # approximation with error bound
    if ɛ < Inf
        M = sparse([1, 2], [block[1], block[2]], [one(N), one(N)], 2, n)
        return overapproximate(M * S, ɛ)
    end

    # box approximation
    plus_one = [one(N)]
    minus_one = [-one(N)]
    pe = σ(sparsevec([block[1]], plus_one, n), S)
    pw = σ(sparsevec([block[1]], minus_one, n), S)
    pn = σ(sparsevec([block[2]], plus_one, n), S)
    ps = σ(sparsevec([block[2]], minus_one, n), S)
    pe_bi = dot(DIR_EAST(N), view(pe, block))
    pn_bi = dot(DIR_NORTH(N), view(pn, block))
    pw_bi = dot(DIR_WEST(N), view(pw, block))
    ps_bi = dot(DIR_SOUTH(N), view(ps, block))
    return HPolygon([LinearConstraint(DIR_EAST(N), pe_bi),
                     LinearConstraint(DIR_NORTH(N), pn_bi),
                     LinearConstraint(DIR_WEST(N), pw_bi),
                     LinearConstraint(DIR_SOUTH(N), ps_bi)])
end

"""
    project(S::LazySet{N},
            block::AbstractVector{Int},
            set_type::Type{<:Hyperrectangle},
            [n]::Int=dim(S),
            [ɛ]::Real=Inf
           )::Hyperrectangle where {N<:Real}

Project a high-dimensional set to a low-dimensional hyperrectangle.

### Input

- `S` -- set
- `block` -- block structure - a vector with the dimensions of interest
- `set_type` -- `Hyperrectangle` - used for dispatch
- `n` -- (optional, default: `dim(S)`) ambient dimension of the set `S`
- `ɛ` -- (optional, default: `Inf`) - used for dispatch, ignored

### Output

The box approximation of the projection of `S`.
"""
@inline function project(S::LazySet{N},
                         block::AbstractVector{Int},
                         set_type::Type{<:Hyperrectangle},
                         n::Int=dim(S),
                         ɛ::Real=Inf
                        )::Hyperrectangle where {N<:Real}
    high = Vector{N}(length(block))
    low = similar(high)
    for i in eachindex(block)
        high[i] = σ(sparsevec([block[i]], [one(N)], n), S)[block[i]]
        low[i] = σ(sparsevec([block[i]], [-one(N)], n), S)[block[i]]
    end
    return Hyperrectangle(high=high, low=low)
end
