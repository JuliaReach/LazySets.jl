"""
    decompose(S::LazySet{N}, set_type::Type=HPolygon
             )::CartesianProductArray where {N<:Real}

Compute an overapproximation of the projections of the given set over each
two-dimensional subspace.

### Input

- `S`        -- set
- `set_type` -- (optional, default: `HPolygon`) type of set approximation in 2D

### Output

A `CartesianProductArray` corresponding to the Cartesian product of
two-dimensional sets of type `set_type`.

### Algorithm

For each 2D block a specific `decompose_2D` method is called, dispatched on the
`set_type` argument.
"""
function decompose(S::LazySet{N}, set_type::Type=HPolygon
                  )::CartesianProductArray where {N<:Real}
    n = dim(S)
    b = div(n, 2)
    result = Vector{set_type{N}}(b)
    @inbounds for bi in 1:b
        result[bi] = decompose_2D(S, n, bi, set_type)
    end
    return CartesianProductArray(result)
end

# polygon with box directions
@inline function decompose_2D(S::LazySet{N}, n::Int, bi::Int,
                              set_type::Type{<:HPolygon}
                             )::HPolygon where {N<:Real}
    pe, pn, pw, ps = box_bounds(S, n, bi)
    block = 2*bi-1:2*bi
    pe_bi = dot(DIR_EAST(N), view(pe, block))
    pn_bi = dot(DIR_NORTH(N), view(pn, block))
    pw_bi = dot(DIR_WEST(N), view(pw, block))
    ps_bi = dot(DIR_SOUTH(N), view(ps, block))

    return HPolygon([LinearConstraint(DIR_EAST(N), pe_bi),
                     LinearConstraint(DIR_NORTH(N), pn_bi),
                     LinearConstraint(DIR_WEST(N), pw_bi),
                     LinearConstraint(DIR_SOUTH(N), ps_bi)])
end

# hyperrectangle
@inline function decompose_2D(S::LazySet, n::Int, bi::Int,
                              set_type::Type{<:Hyperrectangle})::Hyperrectangle
    pe, pn, pw, ps = box_bounds(S, n, bi)
    block = 2*bi-1:2*bi
    radius = [(pe[block[1]] - pw[block[1]]) / 2, (pn[block[2]] - ps[block[2]]) / 2]
    center = [pw[block[1]] + radius[1], ps[block[2]] + radius[2]]
    return Hyperrectangle(center, radius)
end

# helper function
@inline function box_bounds(S::LazySet{N}, n::Int, bi::Int) where {N<:Real}
    pe = σ(sparsevec([2*bi-1], [one(N)], n), S)
    pn = σ(sparsevec([2*bi], [one(N)], n), S)
    pw = σ(sparsevec([2*bi-1], [-one(N)], n), S)
    ps = σ(sparsevec([2*bi], [-one(N)], n), S)
    return (pe, pn, pw, ps)
end

"""
    decompose(S::LazySet{N}, ɛi::Vector{<:Real})::CartesianProductArray where {N<:Real}

Compute an overapproximation of the projections of the given set over each
two-dimensional subspace with a certified error bound.

### Input

- `S`  -- set
- `ɛi` -- array with the error bound for each projection (different error bounds
          can be passed for different blocks)

### Output

A `CartesianProductArray` corresponding to the Cartesian product of ``2 × 2``
polygons.

### Algorithm

This algorithm assumes a decomposition into two-dimensional subspaces only,
i.e., partitions of the form ``[2, 2, …, 2]``.
In particular, if `S` is a `CartesianProductArray`, no check is performed to
verify that assumption.

The algorithm proceeds as follows:

1. Project the set `S` into each partition, with `M⋅S`, where M is the identity
   matrix in the block coordinates and zero otherwise.
2. Overapproximate the set with a given error bound, `ɛi[i]`, for
   ``i = 1, …, b``,
3. Return the result as a `CartesianProductArray`.
"""
function decompose(S::LazySet{N}, ɛi::Vector{<:Real})::CartesianProductArray where {N<:Real}
    n = dim(S)
    b = div(n, 2)
    result = Vector{HPolygon{N}}(b)
    @inbounds for i in 1:b
        M = sparse([1, 2], [2*i-1, 2*i], [one(N), one(N)], 2, n)
        result[i] = overapproximate(M * S, ɛi[i])
    end
    return CartesianProductArray(result)
end

"""
    decompose(S::LazySet, ɛ::Real, [set_type]::Type=HPolygon
             )::CartesianProductArray

Compute an overapproximation of the projections of the given set over each
two-dimensional subspace with a certified error bound.

### Input

- `S`        -- set
- `ɛ`        -- error bound
- `set_type` -- (optional, default: `HPolygon`) type of set approximation in 2D

### Output

A `CartesianProductArray` corresponding to the Cartesian product of
two-dimensional sets of type `set_type`.

### Notes

This function is a particular case of `decompose(S, ɛi)`, where the same error
bound for each block is assumed.

The `set_type` argument is ignored if ``ɛ ≠ \\text{Inf}``.
"""
function decompose(S::LazySet, ɛ::Real, set_type::Type=HPolygon
                  )::CartesianProductArray
    if ɛ == Inf
        return decompose(S, set_type)
    else
        return decompose(S, [ɛ for i in 1:div(dim(S), 2)])
    end
end
