"""
  overapproximate(S)

Return an approximation of the given 2D set as a polygon, using box directions.

INPUT:

- ``S`` -- a 2D set defined by its support function

OUTPUT:

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
    overapproximate(S, ɛ)

Return an ɛ-close approximation of the given 2D set (in terms of Hausdorff
distance) as a polygon.

INPUT:

- ``S`` -- a 2D set defined by its support function
- ``ɛ`` -- the error bound

OUTPUT:

A polygon in constraint representation.
"""
function overapproximate(S::LazySet, ɛ::Float64)::HPolygon
    constraints = [LinearConstraint(la.d1, dot(la.d1, la.p1)) for la in approximate(S, ɛ)]
    return HPolygon(constraints)
end

"""
    decompose(X, ɛi)

Compute an overapproximation of the projections of the given set over each
two-dimensional subspace.

INPUT:

- ``X`` -- set represented by support functions
- ``ɛi`` -- array, error bound for each projection

OUTPUT:

A CartesianProductArray corresponding to the cartesian product of 2x2 polygons.

NOTES:

- It assumes blocks of size 2.
- Different error bounds can be passed to different blocks.
- If the type of the given X is CartesianProductArray, we cannot do simply

```
return CartesianProductArray([overapproximate(X.sfarray[i], ɛi[i]) for i in 1:b])
```
because the cartesian product could contain objects in different dimension, a priori.

ALGORITHM:

1. Project the set X into each partition (assuming two-dimensional subspaces, i.e. partitions of the form [2, 2, ..., 2]),
with M*X, where M is the identity matrix in the block coordinates and zero otherwise.
2. Overapproximate the set with a given error bound, ɛi[i], for i = 1,.., b
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

function decompose(X::LazySet, ɛ::Float64)::CartesianProductArray
    n = LazySets.dim(X)
    b = div(n, 2)
    return decompose(X, [ɛ for i in 1:b])
end
