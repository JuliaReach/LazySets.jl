"""
    overapproximate(S::LazySet)::HPolygon

Return an approximation of a given 2D convex set as a box-shaped polygon.

### Input

- `S` -- convex set, assumed to be two-dimensional

### Output

A box-shaped polygon in constraint representation.
"""
function overapproximate(S::LazySet)::HPolygon
    @assert dim(S) == 2

    # evaluate support vector on box directions
    pe = σ(DIR_EAST, S)
    pn = σ(DIR_NORTH, S)
    pw = σ(DIR_WEST, S)
    ps = σ(DIR_SOUTH, S)
    constraints = Vector{LinearConstraint{eltype(pe)}}(4)
    constraints[1] = LinearConstraint(DIR_EAST, dot(pe, DIR_EAST))
    constraints[2] = LinearConstraint(DIR_NORTH, dot(pn, DIR_NORTH))
    constraints[3] = LinearConstraint(DIR_WEST, dot(pw, DIR_WEST))
    constraints[4] = LinearConstraint(DIR_SOUTH, dot(ps, DIR_SOUTH))
    return HPolygon(constraints)
end

"""
    overapproximate(S::LazySet, ɛ::Float64)::HPolygon

Return an ɛ-close approximation of the given 2D set (in terms of Hausdorff
distance) as a polygon.

### Input

- `S` -- convex set, assumed to be two-dimensional
- `ɛ` -- error bound

### Output

A polygon in constraint representation.
"""
function overapproximate(S::LazySet, ɛ::Float64)::HPolygon
    @assert dim(S) == 2

    constraints =
        [LinearConstraint(ci.d1, dot(ci.d1, ci.p1)) for ci in approximate(S, ɛ)]
    return HPolygon(constraints)
end

"""
    decompose(S::LazySet)::CartesianProductArray

Compute an overapproximation of the projections of the given convex set over
each two-dimensional subspace using box directions.

### Input

- `S` -- convex set

### Output

A `CartesianProductArray` corresponding to the Cartesian product of
``2 \\times 2`` box-shaped polygons.
"""
function decompose(S::LazySet)::CartesianProductArray
    n = dim(S)
    b = div(n, 2)
    result = Vector{HPolygon}(b)

    @inbounds for bi in 1:b
        pe_bi = dot(DIR_EAST,
                    view(σ(sparsevec([2*bi-1], [1.], n), S), 2*bi-1:2*bi))

        pn_bi = dot(DIR_NORTH,
                    view(σ(sparsevec([2*bi], [1.], n), S), 2*bi-1:2*bi))

        pw_bi = dot(DIR_WEST,
                    view(σ(sparsevec([2*bi-1], [-1.], n), S), 2*bi-1:2*bi))

        ps_bi = dot(DIR_SOUTH,
                    view(σ(sparsevec([2*bi], [-1.], n), S), 2*bi-1:2*bi))

        result[bi] = HPolygon([LinearConstraint(DIR_EAST, pe_bi),
                               LinearConstraint(DIR_NORTH, pn_bi),
                               LinearConstraint(DIR_WEST, pw_bi),
                               LinearConstraint(DIR_SOUTH, ps_bi)])
    end

    return CartesianProductArray(result)
end

"""
    decompose(S::LazySet, ɛi::Vector{Float64})::CartesianProductArray

Compute an overapproximation of the projections of the given convex set over
each two-dimensional subspace with a certified error bound.

### Input

- `S`  -- convex set
- `ɛi` -- array with the error bound for each projection (different error bounds
          can be passed for different blocks)

### Output

A `CartesianProductArray` corresponding to the Cartesian product of
``2 \\times 2`` polygons.

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
function decompose(S::LazySet, ɛi::Vector{Float64})::CartesianProductArray
    n = dim(S)
    b = div(n, 2)
    result = Vector{HPolygon}(b)
    @inbounds for i in 1:b
        M = sparse([1, 2], [2*i-1, 2*i], [1., 1.], 2, n)
        result[i] = overapproximate(M * S, ɛi[i])
    end
    return CartesianProductArray(result)
end

"""
    decompose(S::LazySet, ɛ::Float64)::CartesianProductArray

Compute an overapproximation of the projections of the given convex set over
each two-dimensional subspace with a certified error bound.

### Input

- `S` -- convex set
- `ɛ` --  error bound

### Output

A `CartesianProductArray` corresponding to the Cartesian product of
``2 \\times 2`` polygons.

### Notes

This function is a particular case of `decompose(S, ɛi)`, where the same error
bound for each block is assumed.
"""
function decompose(S::LazySet, ɛ::Float64)::CartesianProductArray
    if ɛ == Inf
        return decompose(S)
    else
        return decompose(S, [ɛ for i in 1:div(dim(S), 2)])
    end
end
