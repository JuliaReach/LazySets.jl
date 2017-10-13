"""
  overapproximate(S)

Return an approximation of the given 2D set as a polygon, using box directions.

### Input

- `S` -- a 2D set defined by its support function

### Output

A polygon in constraint representation.
"""
function overapproximate(S::LazySet)::HPolygon
    constraints = Array{LinearConstraint, 1}(4)
    # evaluate support vector on box directions
    pe = σ(DIR_EAST, S)
    pn = σ(DIR_NORTH, S)
    pw = σ(DIR_WEST, S)
    ps = σ(DIR_SOUTH, S)
    constraints[1] = LinearConstraint(DIR_EAST, dot(pe, DIR_EAST))
    constraints[2] = LinearConstraint(DIR_NORTH, dot(pn, DIR_NORTH))
    constraints[3] = LinearConstraint(DIR_WEST, dot(pw, DIR_WEST))
    constraints[4] = LinearConstraint(DIR_SOUTH, dot(ps, DIR_SOUTH))
    return HPolygon(constraints)
end

"""
    overapproximate(S, ɛ)

Return an ɛ-close approximation of the given 2D set (in terms of Hausdorff
distance) as a polygon.

### Input

- `S` -- a 2D set defined by its support function
- `ɛ` -- the error bound

### Output

A polygon in constraint representation.
"""
function overapproximate(X::LazySet, ɛ::Float64)::HPolygon
    constraints = [LinearConstraint(la.d1, dot(la.d1, la.p1)) for la in approximate(X, ɛ)]
    return HPolygon(constraints)
end

"""
    decompose(X)

Compute an overapproximation of the projections of the given set over each
two-dimensional subspace using box directions.

### Input

- `X`  -- set represented by support functions

### Output

A CartesianProductArray corresponding to the cartesian product of 2x2 polygons.
"""
function decompose(X::LazySet)::CartesianProductArray
    n = LazySets.dim(X)
    b = div(n, 2)

    DIR_EAST_b = repeat(DIR_EAST, outer=b)
    DIR_NORTH_b = repeat(DIR_NORTH, outer=b)
    DIR_WEST_b = repeat(DIR_WEST, outer=b)
    DIR_SOUTH_b = repeat(DIR_SOUTH, outer=b)

    pe = reshape(σ(DIR_EAST_b, X), 2, :)
    pe = DIR_EAST.' * pe

    pn = reshape(σ(DIR_NORTH_b, X), 2, :)
    pn = DIR_NORTH.' * pn

    pw = reshape(σ(DIR_WEST_b, X), 2, :)
    pw = DIR_WEST.' * pw

    ps = reshape(σ(DIR_SOUTH_b, X), 2, :)
    ps = DIR_SOUTH.' * ps

    result = Array{HPolygon, 1}(b)
    @inbounds for i in 1:b
        result[i] = HPolygon([LinearConstraint(DIR_EAST, pe[i]),
                              LinearConstraint(DIR_NORTH, pn[i]),
                              LinearConstraint(DIR_WEST, pw[i]),
                              LinearConstraint(DIR_SOUTH, ps[i])])
    end
    return CartesianProductArray(result)
end

"""
    decompose(X, ɛi)

Compute an overapproximation of the projections of the given set over each
two-dimensional subspace with a certified error bound.

### Input

- `X`  -- set represented by support functions
- `ɛi` -- array, error bound for each projection (different error bounds
          can be passed to different blocks)

### Output

A CartesianProductArray corresponding to the cartesian product of 2x2 polygons.

### Algorithm

This algorithm assumes a decomposition into two-dimensional subspaces only,
i.e. partitions of the form ``[2, 2, ..., 2]``. In particular if `X` is a `CartesianProductArray`
no check is performed to verify that assumption.

It proceeds as follows:

1. Project the set `X` into each partition, with ``MX``, where ``M`` is the
identity matrix in the block coordinates and zero otherwise.
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

- `X`  -- set represented by support functions
- `ɛ` --  error bound

### Output

A CartesianProductArray corresponding to the cartesian product of 2x2 polygons.
"""
function decompose(X::LazySet, ɛ::Float64)::CartesianProductArray
    n = LazySets.dim(X)
    b = div(n, 2)
    return decompose(X, [ɛ for i in 1:b])
end
