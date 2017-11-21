"""
  overapproximate(X)

Return an approximation of the given 2D set as a box-shaped polygon.

### Input

- `X` -- lazy set, assumed to be two-dimensional

### Output

A polygon in constraint representation.
"""
function overapproximate(X::LazySet)::HPolygon
    constraints = Vector{LinearConstraint}(4)
    # evaluate support vector on box directions
    pe = σ(DIR_EAST, X)
    pn = σ(DIR_NORTH, X)
    pw = σ(DIR_WEST, X)
    ps = σ(DIR_SOUTH, X)
    constraints[1] = LinearConstraint(DIR_EAST, dot(pe, DIR_EAST))
    constraints[2] = LinearConstraint(DIR_NORTH, dot(pn, DIR_NORTH))
    constraints[3] = LinearConstraint(DIR_WEST, dot(pw, DIR_WEST))
    constraints[4] = LinearConstraint(DIR_SOUTH, dot(ps, DIR_SOUTH))
    return HPolygon(constraints)
end

"""
    overapproximate(X, ɛ)

Return an ɛ-close approximation of the given 2D set (in terms of Hausdorff
distance) as a polygon.

### Input

- `X` -- lazy set, assumed to be two-dimensional
- `ɛ` -- the error bound

### Output

A polygon in constraint representation.
"""
function overapproximate(X::LazySet, ɛ::Float64)::HPolygon
    constraints = [LinearConstraint(ci.d1, dot(ci.d1, ci.p1)) for ci in approximate(X, ɛ)]
    return HPolygon(constraints)
end

"""
    decompose(X)

Compute an overapproximation of the projections of the given set over each
two-dimensional subspace using box directions.

### Input

- `X`  -- lazy set

### Output

A `CartesianProductArray` corresponding to the cartesian product of 2x2 polygons.
"""
function decompose(X::LazySet)::CartesianProductArray
    n = LazySets.dim(X)
    b = div(n, 2)

    DIR_EAST_bi = Vector{Float64}(n)
    DIR_NORTH_bi = Vector{Float64}(n)
    DIR_WEST_bi = Vector{Float64}(n)
    DIR_SOUTH_bi = Vector{Float64}(n)
    result = Vector{HPolygon}(b)

    @inbounds for bi in 1:b
        copy!(DIR_EAST_bi, sparsevec([2*bi-1], [1.], n))
        pe_bi = dot(DIR_EAST, view(σ(DIR_EAST_bi, X), 2*bi-1:2*bi))

        copy!(DIR_NORTH_bi, sparsevec([2*bi], [1.], n))
        pn_bi = dot(DIR_NORTH, view(σ(DIR_NORTH_bi, X), 2*bi-1:2*bi))

        copy!(DIR_WEST_bi, sparsevec([2*bi-1], [-1.], n))
        pw_bi = dot(DIR_WEST, view(σ(DIR_WEST_bi, X), 2*bi-1:2*bi))

        copy!(DIR_SOUTH_bi, sparsevec([2*bi], [-1.], n))
        ps_bi = dot(DIR_SOUTH, view(σ(DIR_SOUTH_bi, X), 2*bi-1:2*bi))

        result[bi] = HPolygon([LinearConstraint(DIR_EAST, pe_bi),
                               LinearConstraint(DIR_NORTH, pn_bi),
                               LinearConstraint(DIR_WEST, pw_bi),
                               LinearConstraint(DIR_SOUTH, ps_bi)])
    end

    return CartesianProductArray(result)
end

"""
    decompose(X, ɛi)

Compute an overapproximation of the projections of the given set over each
two-dimensional subspace with a certified error bound.

### Input

- `X`  -- lazy set
- `ɛi` -- array with the error bound for each projection (different error bounds
          can be passed to different blocks)

### Output

A `CartesianProductArray` corresponding to the cartesian product of 2x2 polygons.

### Algorithm

This algorithm assumes a decomposition into two-dimensional subspaces only,
i.e. partitions of the form ``[2, 2, ..., 2]``.
In particular if `X` is a `CartesianProductArray`, no check is performed
to verify that assumption.

It proceeds as follows:

1. Project the set `X` into each partition, with `MX`, where M is the identity
   matrix in the block coordinates and zero otherwise.
2. Overapproximate the set with a given error bound, `ɛi[i]`, for ``i = 1, …, b``,
3. Return the result as an array of support functions.
"""
function decompose(X::LazySet, ɛi::Vector{Float64})::CartesianProductArray
    n = LazySets.dim(X)
    b = div(n, 2)
    result = Array{HPolygon, 1}(b)
    @inbounds for i in 1:b
        M = sparse([1, 2], [2*i-1, 2*i], [1., 1.], 2, n)
        result[i] = overapproximate(M * X, ɛi[i])
    end
    return CartesianProductArray(result)
end

"""
    decompose(X, ɛ)

Compute an overapproximation of the projections of the given set over each
two-dimensional subspace with a certified error bound.

This function is a particular case of `decompose(X, ɛi)`, where the same error
bound for each block is assumed.

### Input

- `X`  -- lazy set
- `ɛ` --  error bound

### Output

A `CartesianProductArray` corresponding to the cartesian product of 2x2 polygons.
"""
function decompose(X::LazySet, ɛ::Float64)::CartesianProductArray
    n = LazySets.dim(X)
    b = div(n, 2)
    return decompose(X, [ɛ for i in 1:b])
end
